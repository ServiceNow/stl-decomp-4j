package com.github.servicenow.ds.stats.stl;

/**
 * LoessSmoother uses LOESS interpolation to compute a smoothed data set from a regularly-spaced set of input
 * data. If a jump is specified, then LOESS interpolation is only done on every jump points and linear interpolation is
 * done to fill in the gaps.
 *
 * Author: Jim Crotinger, ported from the original RATFOR source from netlib
 */
@SuppressWarnings("WeakerAccess")
public class LoessSmoother {

	private final LoessInterpolator fInterpolator;
	private final double[] fData;
	private final int fWidth;
	private final int fJump;
	private final double[] fSmoothed;

	public static class Builder {
		private Integer fWidth = null;
		private int fDegree = 1;
		private int fJump = 1;
		private double[] fExternalWeights = null;
		private double[] fData = null;

		public Builder setWidth(int width) {
			fWidth = width;
			return this;
		}

		public Builder setDegree(int degree) {
			if (degree < 0 || degree > 2)
				throw new IllegalArgumentException("Degree must be 0, 1 or 2");

			fDegree = degree;
			return this;
		}

		public Builder setJump(int jump) {
			fJump = jump;
			return this;
		}

		public Builder setExternalWeights(double[] weights) {
			fExternalWeights = weights;
			return this;
		}

		public Builder setData(double[] data) {
			fData = data;
			return this;
		}

		public LoessSmoother build() {
			if (fWidth == null)
				throw new IllegalStateException("LoessSmoother.Builder: Width must be set before calling build");

			if (fData == null)
				throw new IllegalStateException("LoessSmoother.Builder: Data must be set before calling build");

			return new LoessSmoother(fWidth, fJump, fDegree, fData, fExternalWeights);
		}
	}

	// -----------------------------------------------------------------------------------------------------------------
	// Interface
	// -----------------------------------------------------------------------------------------------------------------

	/**
	 * Create a LoessSmoother for the given data set with the specified smoothing width and optional external
	 * Weights.
	 *
	 * @param width
	 *            int approximate width the width of the neighborhood weighting function
	 * @param jump
	 *            int smoothing jump - only ever jump points are smoothed by LOESS with linear interpolation in between.
	 * @param degree
	 *            int 1 for linear regression, 0 for simple weighted average
	 * @param data
	 *            double[] underlying data set that is being smoothed
	 * @param externalWeights
	 *            double[] additional weights to apply in the smoothing. Ignored if null.
	 */
	private LoessSmoother(int width, int jump, int degree, double[] data, double[] externalWeights) {
		final LoessInterpolator.Builder b = new LoessInterpolator.Builder();
		this.fInterpolator = b.setWidth(width).setDegree(degree).setExternalWeights(externalWeights).interpolate(data);
		this.fData = data;
		this.fJump = Math.min(jump, data.length - 1);
		this.fWidth = width;
		this.fSmoothed = new double[data.length];
	}

	/**
	 * Accessor to retrieve the underlying interpolator.
	 *
	 * @return LoessInterpolator the underlying interpolators
	 */
	public LoessInterpolator getInterpolator() {
		return fInterpolator;
	}

	// TODO: Refactor to use a strategy pattern - dependencies on final params are determined at construction time.

	/**
	 * Calculate the LOESS smoothed data for each original data point.
	 *
	 * @return double[] array containing the results
	 */
	public double[] smooth() {
		if (fData.length == 1) {
			fSmoothed[0] = fData[0];
			return fSmoothed;
		}

		int left = -1, right = -1;
		if (fWidth >= fData.length) {
			left = 0;
			right = fData.length - 1;
			for (int i = 0; i < fData.length; i += fJump) {
				final Double y = fInterpolator.smoothOnePoint(i, left, right);
				fSmoothed[i] = y == null ? fData[i] : y;
				// logSmoothedPoint(i, smooth[i]);
			}
		} else if (fJump == 1) {
			final int halfWidth = (fWidth + 1) / 2;
			left = 0;
			right = fWidth - 1;
			for (int i = 0; i < fData.length; ++i) {
				if (i >= halfWidth && right != fData.length - 1) {
					++left;
					++right;
				}
				final Double y = fInterpolator.smoothOnePoint(i, left, right);
				fSmoothed[i] = y == null ? fData[i] : y;
				// logSmoothedPoint(i, smooth[i]);
			}
		} else {
			// For reference, the original RATFOR:
			// else { # newnj greater than one, len less than n
			// nsh = (len+1)/2
			// do i = 1,n,newnj { # fitted value at i
			// if(i<nsh) {              // i     = [1, 2, 3, 4, 5, 6, 7, 8, 9]; 9 points
			// nleft = 1                // left  = [1, 1, 1, 1, 1, 1, 1, 1, 1];
			// nright = len             // right = [19, 19, 19, 19, 19, 19, 19, 19, 19]; right - left = 18
			// }
			// else if(i>=n-nsh+1) {    // i     = [135, 136, 137, 138, 139, 140, 141, 142, 143, 144]; 10 points
			// nleft = n-len+1          // left  = [126, 126, 126, 126, 126, 126, 126, 126, 126, 126];
			// nright = n               // right = [144, 144, 144, 144, 144, 144, 144, 144, 144, 144]; right - left = 18
			// }
			// else {                   // i     = [10, 11, 12, ..., 132, 133, 134]; 125 points
			// nleft = i-nsh+1          // left  = [1, 2, 3, ..., 123, 124, 125]
			// nright = len+i-nsh       // right = [19, 20, 21, ..., 141, 142, 143]; right - left = 18
			// }
			// call est(y,n,len,ideg,float(i),ys(i),nleft,nright,res,userw,rw,ok)
			// if(!ok) ys(i) = y(i)
			// }
			// }
			// Note that RATFOR/Fortran are indexed from 1
			//
			// test: data.length == 144, fWidth = 19
			//   --> halfWidth = 10
			// Ignoring jumps...
			// First branch for  i = [0, 1, 2, 3, 4, 5, 6, 7, 8]; 9 points
			//                left = [0, 0, 0, 0, 0, 0, 0, 0, 0]
			//               right = [18, 18, 18, 18, 18, 18, 18, 18, 18]; right - left = 18
			// Second branch for i = [134, 135, 136, 137, 138, 139, 140, 141, 142, 143]; 10 points
			//                left = [125, 125, 125, 125, 125, 125, 125, 125, 125, 125];
			//               right = [143, 143, 143, 143, 143, 143, 143, 143, 143, 143]; right - left = 18
			// Third branch for  i = [ 9, 10, 11, ..., 131, 132, 133]; 125 points
			//                left = [ 0,  1,  2, ..., 122, 123, 124]
			//               right = [18, 19, 20, ..., 140, 141, 142]; right - left = 18
			final int halfWidth = (fWidth + 1) / 2;
			for (int i = 0; i < fData.length; i += fJump) {
				if (i < halfWidth - 1) {
					left = 0;
				} else if (i >= fData.length - halfWidth) {
					left = fData.length - fWidth;
				} else {
					left = i - halfWidth + 1;
				}
				right = left + fWidth - 1;
				final Double y = fInterpolator.smoothOnePoint(i, left, right);
				fSmoothed[i] = y == null ? fData[i] : y;
				// logSmoothedPoint(i, smooth[i]);
			}
		}

		if (fJump != 1) {
			for (int i = 0; i < fData.length - fJump; i += fJump) {
				final double slope = (fSmoothed[i + fJump] - fSmoothed[i]) / (double) fJump;
				for (int j = i + 1; j < i + fJump; ++j) {
					fSmoothed[j] = fSmoothed[i] + slope * (j - i);
					// logInterpolatedPoint(j, smooth[j]);
				}
			}

			final int last = fData.length - 1;
			int lastSmoothedPos = (last / fJump) * fJump;
			if (lastSmoothedPos != last) {
				final Double y = fInterpolator.smoothOnePoint(last, left, right);
				fSmoothed[last] = y == null ? fData[last] : y;
				// logSmoothedPoint(last, smooth[last]);

				if (lastSmoothedPos != last - 1) {
					final double slope = (fSmoothed[last] - fSmoothed[lastSmoothedPos]) / (last - lastSmoothedPos);
					for (int j = lastSmoothedPos + 1; j < last; ++j) {
						fSmoothed[j] = fSmoothed[lastSmoothedPos] + slope * (j - lastSmoothedPos);
						// logInterpolatedPoint(j, smooth[j]);
					}
				}
			}
		}

		return fSmoothed;
	}

	// private void logSmoothedPoint(double x, double y) {
	// System.out.println(String.format("Smoothed value y(%f) = %f", x, y));
	// }
	//
	// private void logInterpolatedPoint(double x, double y) {
	// System.out.println(String.format("Linear interpolated value y(%f) = %f", x, y));
	// }
}
