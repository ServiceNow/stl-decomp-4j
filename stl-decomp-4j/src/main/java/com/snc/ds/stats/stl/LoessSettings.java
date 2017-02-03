package com.snc.ds.stats.stl;

/**
 * LoessSettings - immutable settings class for specifying the triple of width/jump/degree used to initialize the LOESS
 * interpolators. Enforces that width be 3 or larger and that the degree is 0 or 1.
 *
 * Created by Jim Crotinger on 12-May-2016.
 *
 */
@SuppressWarnings("WeakerAccess")
public class LoessSettings {

	private final int fWidth;
	private final int fDegree;
	private final int fJump;

	/**
	 * Create specifying width, degree and jump.
	 *
	 * @param width
	 *            int width of the LOESS smoother in data points.
	 * @param degree
	 *            int degree of polynomial used in LOESS.
	 * @param jump
	 *            int number of points to skip between LOESS smoothings.
	 */
	public LoessSettings(int width, int degree, int jump) {
		width = Math.max(3, width);
		if (width % 2 == 0)
			++width;
		this.fWidth = width;
		this.fJump = Math.max(1, jump);
		degree = Math.max(0, Math.min(2, degree));
		this.fDegree = degree;
	}

	/**
	 * Create specifying width and degree, defaulting jump to 10% of smoothing width.
	 *
	 * @param width
	 *            int width of the LOESS smoother in data points.
	 * @param degree
	 *            int degree of polynomial used in LOESS.
	 */
	public LoessSettings(int width, int degree) {
		// NOTE: calling this(width, degree, Math.max(1, (int) (0.1 * width + 0.9))) is wrong here since width hasn't
		// been adjusted yet. Simpler to just copy the code and test.
		width = Math.max(3, width);
		if (width % 2 == 0)
			++width;
		this.fWidth = width;
		this.fJump = Math.max(1, (int) (0.1 * width + 0.9));
		degree = Math.max(0, Math.min(2, degree));
		this.fDegree = degree;
	}

	/**
	 * Create from width only. Defaults to linear degree and a jump of 10% of the smoothing width. Enforces minimal
	 * width of 3.
	 *
	 * @param width
	 *            int width of the LOESS smoother in data points.
	 */
	public LoessSettings(int width) {
		width = Math.max(3, width);
		if (width % 2 == 0)
			++width;
		this.fWidth = width;
		this.fJump = Math.max(1, (int) (0.1 * width + 0.9));
		this.fDegree = 1;
	}

	/**
	 * Get the width of the LOESS smoother.
	 *
	 * @return int width of LOESS smoother.
	 */
	public final int getWidth() {
		return fWidth;
	}

	/**
	 * Get the degree of the LOESS smoother
	 *
	 * @return int degree of the LOESS smoother.
	 */
	public final int getDegree() {
		return fDegree;
	}

	/**
	 * Get the jump used between LOESS interpolations.
	 *
	 * @return int jump used between LOESS interpolations.
	 */
	public final int getJump() {
		return fJump;
	}

	@Override
	public String toString() {
		return String.format("[width = %d, degree = %d, jump = %d]", fWidth, fDegree, fJump);
	}
}
