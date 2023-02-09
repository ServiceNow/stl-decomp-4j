package com.github.servicenow.ds.stats.stl;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

/**
 * LoessInterpolator implements Cleveland et al's LOESS smoothing from their STL functionality.
 * <p>
 * Author: Jim Crotinger, based on the original RATFOR source from netlib, with quadratic regression added May 2016.
 */
abstract class LoessInterpolator {

	// -----------------------------------------------------------------------------------------------------------------
	// State
	// -----------------------------------------------------------------------------------------------------------------

	private final int fWidth;
	private final double[] fExternalWeights;

	double[] fData;
	double[][] fExogenousData;
	final double[] fWeights;
	final boolean fOutputNonExogenousPart;

	// -----------------------------------------------------------------------------------------------------------------
	// Construction
	// -----------------------------------------------------------------------------------------------------------------

	/**
	 * class LoessInterpolator.Builder - Factory for LoessInterpolator objects
	 * <p>
	 * Use a Builder pattern to support default arguments and eliminate the multiplicity of same-typed arguments.
	 */
	static class Builder {

		private Integer fWidth = null;
		private int fDegree = 1;
		private double[] fExternalWeights = null;
		private boolean fOutputNonExogenousPart = false;

		/**
		 * Set the width of the LOESS smoother.
		 *
		 * @param width
		 * @return this
		 */
		public Builder setWidth(int width) {
			fWidth = width;
			return this;
		}

		/**
		 * Set the degree of the LOESS interpolator. Defaults to 1.
		 *
		 * @param degree
		 * @return this
		 */
		Builder setDegree(int degree) {
			if (degree < 0 || degree > 2)
				throw new IllegalArgumentException("Degree must be 0, 1 or 2");

			fDegree = degree;
			return this;
		}

		/**
		 * Set the external weights for interpolation.
		 *
		 * Not required - null is equivalent to all weights being 1.
		 *
		 * @param weights
		 * @return this
		 */
		Builder setExternalWeights(double[] weights) {
			fExternalWeights = weights;
			return this;
		}

		/**
		 * Set the boolean for outputting only the nonexogenous part from the smoother.
		 *
		 * @param outputNonExogenousPart boolean if True then only outputs the const+trend from the smoother.
		 * @return this
		 */
		public Builder setOutputNonExogenousPart(boolean outputNonExogenousPart) {
			fOutputNonExogenousPart = outputNonExogenousPart;
			return this;
		}

		/**
		 * Create a LoessInterpolator for interpolating the given data array.
		 *
		 * @param data
		 * @param exogenousData
		 * @return new LoessInterpolator
		 */
		public LoessInterpolator interpolate(double[] data, double[][] exogenousData) {
			if (fWidth == null)
				throw new IllegalStateException("LoessInterpolator.Builder: Width must be set");

			if (data == null)
				throw new IllegalStateException("LoessInterpolator.Builder: data must be non-null");

			switch (fDegree) {
				case 0:
					return new FlatLoessInterpolator(fWidth, data, exogenousData, fExternalWeights, fOutputNonExogenousPart);
				case 1:
					return new LinearLoessInterpolator(fWidth, data, exogenousData, fExternalWeights, fOutputNonExogenousPart);
				case 2:
					return new QuadraticLoessInterpolator(fWidth, data, exogenousData, fExternalWeights, fOutputNonExogenousPart);
				default:
					return null; // Can't actually get here but compiler didn't figure that out.
			}
		}

		LoessInterpolator interpolate(double[] data) {
			return interpolate(data, null);
		}

	}

	// -----------------------------------------------------------------------------------------------------------------
	// Interface
	// -----------------------------------------------------------------------------------------------------------------

	/**
	 * Given a set of data on the regular grid {left, left+1, ..., right-1, right}, computed the LOESS-smoothed value at
	 * the position x and return it. If the value can't be computed, return null.
	 *
	 * @param x
	 *            double x-coordinate at which we want to compute an estimate of y
	 * @param left
	 *            int leftmost coordinate to use from the input data
	 * @param right
	 *            int rightmost coordinate to use from the input data
	 * @return Double interpolated value, or null if interpolation could not be done
	 */
	Double smoothOnePoint(final double x, final int left, final int right) {

		// Ordinarily, one doesn't do linear regression one x-value at a time, but LOESS does since
		// each x-value will typically have a different window. As a result, the weighted linear regression
		// is recast as a linear operation on the input data, weighted by this.fWeights.

		State state = computeNeighborhoodWeights(x, left, right);

		if (state == State.WEIGHTS_FAILED)
			return null;

		if (fExogenousData != null && state == State.LINEAR_OK) {
			return smoothOnePointExogModel(x, left, right, 1, fOutputNonExogenousPart);
		}

		if (state == State.LINEAR_OK)
			updateWeights(x, left, right);

		double ys = 0.0;
		for (int i = left; i <= right; ++i)
			ys += fWeights[i] * fData[i];

		return ys;
	}

	/**
	 * Update the weights for the appropriate least-squares interpolation.
	 *
	 * @param x
	 *            double x-coordinate at which we want to compute an estimate of y
	 * @param left
	 *            int leftmost coordinate to use from the input data
	 * @param right
	 *            int rightmost coordinate to use from the input data
	 */
	abstract void updateWeights(double x, int left, int right);


	// -----------------------------------------------------------------------------------------------------------------
	// Implementation
	// -----------------------------------------------------------------------------------------------------------------

	/**
	 * Internal enum used to return the state of the weights calculation.
	 */
	private enum State {
		WEIGHTS_FAILED, LINEAR_FAILED, LINEAR_OK
	}

	/**
	 * Computer the neighborhood weights.
	 *
	 * @param x
	 *            double x-coordinate at which we want to compute an estimate of y
	 * @param left
	 *            int leftmost coordinate to use from the input data
	 * @param right
	 *            int rightmost coordinate to use from the input data
	 * @return State indicating the whether we can do linear, moving average, or nothing
	 */
	private State computeNeighborhoodWeights(double x, int left, int right) {

		double lambda = Math.max(x - left, right - x);

		// Ordinarily, lambda ~ width / 2.
		//
		// If width > n, then we will only be computing with n points (i.e. left and right will always be in the
		// domain of 1..n) and the above calculation will give lambda ~ n / 2. We want the shape of the neighborhood
		// weight function to be driven by width, not by the size of the domain, so we adjust lambda to be ~ width / 2.
		// (The paper does this by multiplying the above lambda by (width / n). Not sure why the code is different.)

		if (fWidth > fData.length)
			lambda += (double) ((fWidth - fData.length) / 2);

		// "Neighborhood" is computed somewhat fuzzily.

		final double l999 = 0.999 * lambda;
		final double l001 = 0.001 * lambda;

		// Compute neighborhood weights, updating with external weights if supplied.

		double totalWeight = 0.0;
		for (int j = left; j <= right; ++j) {
			final double delta = Math.abs(x - j);

			// Compute the tri-cube neighborhood weight

			double weight = 0.0;
			if (delta <= l999) {
				if (delta <= l001) {
					weight = 1.0;
				} else {
					final double fraction = delta / lambda;
					final double trix = 1.0 - fraction * fraction * fraction;
					weight = trix * trix * trix;
				}

				// If external weights are provided, apply them.

				if (fExternalWeights != null)
					weight *= fExternalWeights[j];

				totalWeight += weight;
			}

			fWeights[j] = weight;
		}

		// If the total weight is 0, we can't proceed, so signal failure.

		if (totalWeight <= 0.0)
			return State.WEIGHTS_FAILED;

		// Normalize the weights

		for (int j = left; j <= right; ++j)
			fWeights[j] /= totalWeight;

		return (lambda > 0) ? State.LINEAR_OK : State.LINEAR_FAILED;
	}

	/**
	 * Create a LoessInterpolator interpolator for the given data set with the specified smoothing width and
	 * optional external Weights.
	 *
	 * @param width
	 *            - the width of the neighborhood weighting function
	 * @param data
	 *            - underlying data set that is being smoothed
	 * @param exogenousData
	 *            -
	 * @param externalWeights
	 *            - additional weights to apply in the smoothing. Ignored if null.
	 * @param outputNonExogenousPart
	 *            -
	 */
	LoessInterpolator(int width, double[] data,double[][] exogenousData, double[] externalWeights, boolean outputNonExogenousPart) {
		this.fWidth = width;
		this.fData = data;
		this.fExogenousData = exogenousData;
		this.fExternalWeights = externalWeights;
		this.fOutputNonExogenousPart = outputNonExogenousPart;
		this.fWeights = new double[data.length];
	}

	public void setExogenousInputs(double[][] exogenousData) {
		fExogenousData = exogenousData;
	}

	public void setData(double[] data) {
		fData = data;
	}

	protected final double smoothOnePointExogModel(double x, int left, int right, int degree, boolean bOutputNonExogenousPart) {
		// Update weights routine when there are exogenous inputs (solved via weighed least-squares)

		// First form the regressor matrix and the data vector
		int datalength = fData.length;
		int windowlength = right - left + 1;
		int numExogInputs = fExogenousData.length;
		int numRegressorColumns = degree + 1 + numExogInputs;
		double[] xpoints = IntStream.rangeClosed(left, right).mapToDouble(el->el/(double) datalength).toArray();

		double[][] regressorMatrixTranspose = new double[numRegressorColumns][windowlength];
		regressorMatrixTranspose[0] = DoubleStream.generate(()->1.0).limit(windowlength).toArray();
		for (int i = 1; i <=degree; ++i) {
			int intpower = i;
			regressorMatrixTranspose[i] = DoubleStream.of(xpoints).map(el->Math.pow(el, intpower)).toArray();
		}
		for (int i = degree+1; i <numRegressorColumns; ++i) {
			regressorMatrixTranspose[i] = Arrays.copyOfRange(fExogenousData[i-degree-1], left, right+1);
		}
		double[][] regressorMatrix = MatrixUtils.createRealMatrix(regressorMatrixTranspose).transpose().getData();
		double[] datavector = Arrays.copyOfRange(fData, left, right+1);
		double[] weightsvector = Arrays.copyOfRange(fWeights, left, right+1);

		// Then solve the weighted least-squares problem
		double[] datavectorweighted = new double[windowlength];
		double[][] regressorMatrixWeighted = new double[windowlength][numRegressorColumns];

		for (int i = 0; i < windowlength; ++i) {
			double weight = Math.sqrt( Math.max(Math.abs(weightsvector[i]), 1e-20) );
			datavectorweighted[i] = weight * datavector[i];
			for (int j = 0; j < numRegressorColumns; ++j)
				regressorMatrixWeighted[i][j] = weight * regressorMatrix[i][j];
		}
		double[] parameters = leastSquaresEstimation(regressorMatrixWeighted, datavectorweighted);

		int indexRegressors = bOutputNonExogenousPart ? degree + 1 : numRegressorColumns;
		double ys = 0.0;
		for (int i = 0; i < indexRegressors; ++i)
			ys += regressorMatrix[(int) (x-left)][i] * parameters[i];
		return ys;

	}

	protected final double[] leastSquaresEstimation(double[][] regressorMatrix, double[] datavector) {
		/**
		 * Ordinary least-squares (first tries normal inverse, if errors appear, then switches to the pseudoinverse via SVD
		 * @param regressorMatrix the regressor matrix
		 * @param datavector the datavector that is fit
		 * @return vector of estimated parameters
		 */
		try {
			OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
			regression.setNoIntercept(true);
			regression.newSampleData(datavector, regressorMatrix);
			return regression.estimateRegressionParameters();
		}
		// use the pseudoinverse in case of errors
		catch (Exception e) {
			SingularValueDecomposition svd = new SingularValueDecomposition(MatrixUtils.createRealMatrix(regressorMatrix));
			DecompositionSolver solver = svd.getSolver();
			return solver.solve(MatrixUtils.createRealVector(datavector)).toArray();
		}
	}
}

class FlatLoessInterpolator extends LoessInterpolator {
	/**
	 * Create a LoessInterpolator interpolator for the given data set with the specified smoothing width and
	 * optional external Weights.
	 *
	 * @param width           - the width of the neighborhood weighting function
	 * @param data            - underlying data set that is being smoothed
	 * @param externalWeights
	 */
	FlatLoessInterpolator(int width, double[] data, double[][] exogenousData, double[] externalWeights, boolean outputNonExogenousPart) {
		super(width, data, exogenousData, externalWeights, outputNonExogenousPart);
	}

	/**
	 * Weight update for FlatLinearInterpolator is a no-op.
	 */
	final protected void updateWeights(double x, int left, int right) {
	}
}

class LinearLoessInterpolator extends LoessInterpolator {

	/**
	 * Create a LoessInterpolator interpolator for the given data set with the specified smoothing width and
	 * optional external Weights.
	 *
	 * @param width           - the width of the neighborhood weighting function
	 * @param data            - underlying data set that is being smoothed
	 * @param externalWeights
	 */
	LinearLoessInterpolator(int width, double[] data, double[][] exogenousData, double[] externalWeights, boolean outputNonExogenousPart) {
		super(width, data,exogenousData, externalWeights, outputNonExogenousPart);
	}

	/**
	 * Compute weighted least squares fit to the data points and adjust the weights with the results.
	 *
	 * @param x
	 *            double x-coordinate at which we want to compute an estimate of y
	 * @param left
	 *            int leftmost coordinate to use from the input data
	 * @param right
	 *            int rightmost coordinate to use from the input data
	 */
	final protected void updateWeights(double x, int left, int right) {
		double xMean = 0.0;
		for (int i = left; i <= right; ++i)
			xMean += i * fWeights[i];

		double x2Mean = 0.0;
		for (int i = left; i <= right; ++i) {
			final double delta = i - xMean;
			x2Mean += fWeights[i] * delta * delta;
		}

		// Finding y(x) from the least-squares fit can be cast as a linear operation on the input data.
		// This is implemented by updating the weights to include the least-squares weighting of the points.
		// Note that this is only done if the points are spread out enough (variance > (0.001 * range)^2)
		// to compute a slope. If not, we leave the weights alone and essentially fall back to a moving
		// average of the data based on the neighborhood and external weights.

		final double range = fData.length - 1;
		if (x2Mean > 0.000001 * range * range) {
			final double beta = (x - xMean) / x2Mean;
			for (int i = left; i <= right; ++i)
				fWeights[i] *= (1.0 + beta * (i - xMean));
		}
	}
}

class QuadraticLoessInterpolator extends LoessInterpolator {

	/**
	 * Create a QuadraticLoessInterpolator interpolator for the given data set with the specified smoothing width and
	 * optional external Weights.
	 *
	 * @param width           - the width of the neighborhood weighting function
	 * @param data            - underlying data set that is being smoothed
	 * @param externalWeights
	 */
	QuadraticLoessInterpolator(int width, double[] data, double[][] exogenousData, double[] externalWeights, boolean outputNonExogenousPart) {
		super(width, data, exogenousData, externalWeights, outputNonExogenousPart);
	}

	/**
	 * Compute weighted least squares quadratic fit to the data points and adjust the weights with the results.
	 *
	 * @param x
	 *            double x-coordinate at which we want to compute an estimate of y
	 * @param left
	 *            int leftmost coordinate to use from the input data
	 * @param right
	 *            int rightmost coordinate to use from the input data
	 */
	final protected void updateWeights(double x, int left, int right) {

		double x1Mean = 0.0;
		double x2Mean = 0.0;
		double x3Mean = 0.0;
		double x4Mean = 0.0;
		for (int i = left; i <= right; ++i) {
			final double w = fWeights[i];
			final double x1w = i * w;
			final double x2w = i * x1w;
			final double x3w = i * x2w;
			final double x4w = i * x3w;
			x1Mean += x1w;
			x2Mean += x2w;
			x3Mean += x3w;
			x4Mean += x4w;
		}

		final double m2 = x2Mean - x1Mean * x1Mean;
		final double m3 = x3Mean - x2Mean * x1Mean;
		final double m4 = x4Mean - x2Mean * x2Mean;

		final double denominator = m2 * m4 - m3 * m3;
		final double range = fData.length - 1;

		if (denominator > 0.000001 * range * range) {
			// TODO: Are there cases where denominator is too small but m2 is not too small?
			// In that case, it would make sense to fall back to linear regression instead of falling back to just the
			// weighted average.
			final double beta2 = m4 / denominator;
			final double beta3 = m3 / denominator;
			final double beta4 = m2 / denominator;

			final double x1 = x - x1Mean;
			final double x2 = x*x - x2Mean;

			final double a1 = beta2 * x1 - beta3 * x2;
			final double a2 = beta4 * x2 - beta3 * x1;

			for (int i = left; i <= right; ++i)
				fWeights[i] *= (1 + a1 * (i - x1Mean) + a2 * (i * i - x2Mean));

		}
	}
}