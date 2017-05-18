package com.github.servicenow.ds.stats.stl;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * While this class has no behavior, per se, the constructors are robust and fix bad args. Make sure they're doing it
 * right and doing it consistently.
 */
public class LoessSettingsTest {

	@Test
	public void evenWidthBecomesNextOdd() {
		LoessSettings settings = new LoessSettings(20);
		assertEquals(21, settings.getWidth());
		assertEquals(1, settings.getDegree());
		assertEquals(3, settings.getJump());
	}

	/**
	 * Test that all constructors work consistently
	 */
	@Test
	public void evenWidthBecomesNextOdd2() {
		LoessSettings settings = new LoessSettings(20, 0);
		assertEquals(21, settings.getWidth());
		assertEquals(0, settings.getDegree());
		assertEquals(3, settings.getJump());
	}

	/**
	 * Test that all constructors work consistently
	 */
	@Test
	public void evenWidthBecomesNextOdd3() {
		LoessSettings settings = new LoessSettings(20, 0, 4);
		assertEquals(21, settings.getWidth());
		assertEquals(0, settings.getDegree());
		assertEquals(4, settings.getJump());
	}

	@Test
	public void defaultJumpCalculationIsConsistentForOddWidth() {
		LoessSettings settings1 = new LoessSettings(51, 0);
		LoessSettings settings2 = new LoessSettings(51);
		assertEquals(6, settings1.getJump());
		assertEquals(6, settings2.getJump());
	}

	/**
	 * Test for bug where jump was calculated before width was made odd.
	 */
	@Test
	public void defaultJumpCalculationIsConsistentForEvenWidth() {
		LoessSettings settings1 = new LoessSettings(50, 0);
		LoessSettings settings2 = new LoessSettings(50);
		assertEquals(6, settings1.getJump());
		assertEquals(6, settings2.getJump());
	}

	@Test
	public void minWidthIsThree() {
		LoessSettings settings = new LoessSettings(0);
		assertEquals(3, settings.getWidth());
		assertEquals(1, settings.getDegree());
		assertEquals(1, settings.getJump());
	}

	@Test
	public void jumpIsCorrect() {
		LoessSettings settings = new LoessSettings(100);
		assertEquals(11, settings.getJump());
		assertEquals(1, settings.getDegree());
		assertEquals(101, settings.getWidth());
	}

	/**
	 * Test that cap was fixed after LOESS was extended to quadratic - this was broken for some time.
	 */
	@Test
	public void quadraticWorks() {
		LoessSettings settings = new LoessSettings(13, 2, 1);
		assertEquals(2, settings.getDegree());
		assertEquals(13, settings.getWidth());
		assertEquals(1, settings.getJump());
	}

	@Test
	public void jumpIsFlooredAtOne() {
		LoessSettings settings = new LoessSettings(13, 2, -1);
		assertEquals(2, settings.getDegree());
		assertEquals(13, settings.getWidth());
		assertEquals(1, settings.getJump());
	}

	@Test
	public void degreeIsFlooredAtZero() {
		LoessSettings settings = new LoessSettings(13, -2);
		assertEquals(13, settings.getWidth());
		assertEquals(0, settings.getDegree());
		assertEquals(2, settings.getJump());
	}

	@Test
	public void degreeIsCappedAt2() {
		LoessSettings settings = new LoessSettings(13, 10);
		assertEquals(13, settings.getWidth());
		assertEquals(2, settings.getDegree());
		assertEquals(2, settings.getJump());
	}

	@Test
	public void toStringTest() {
		LoessSettings settings = new LoessSettings(23);
		String str = settings.toString();
		assertEquals("[width = 23, degree = 1, jump = 3]", str);
	}
}
