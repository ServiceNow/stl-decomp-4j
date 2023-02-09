package com.snc.ds.stats.stl;

import org.junit.Test;

import static org.junit.Assert.*;

import java.util.Arrays;

/**
 * Regression tests against results generated by original Java implementation.
 */
public class StlMultivariateTest {
    // 2 tests of the multivariate STL with synthetic data

	@Test
	public void RegressionTest1() {
	    // Tests a single exogenous input of a shift occurring at index 70 onwards, 
	    // with indices up to 80 taken as the data, with the periodicity of 4.

        double[] data = Arrays.copyOfRange(fTestData1, 0, 80);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();
		builder.setPeriodLength(4);
		builder.setSeasonalWidth(5);
		builder.setInnerIterations(1);
		builder.setRobustnessIterations(15);

		SeasonalTrendLoess smoother = builder.buildSmoother(data, fExogenousInputs1);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] seasonal = stl.getSeasonal();
		double[] trendexog = stl.getTrend();
		double[] residual = stl.getResidual();
		double epsilon = 1.0e-2;
  
	    assertArrayEquals(fSeas1, seasonal, epsilon);
        assertArrayEquals(fTrendExog1, trendexog, epsilon);
        assertArrayEquals(fResidual1, residual, epsilon);
	}
	
    @Test
    public void RegressionTest2() {
        // Tests two exogenous inputs, first one has instantaneous jumps (outliers) occurring at indices 5, 21, 90
        // whereas the second exogenous input is built from an order 3 autoregressive process that acts as an auxiliary time-series
        // only indices up to 86 are provided as the data, with the periodicity of 7.

        double[] data = Arrays.copyOfRange(fTestData2, 0, 86);


        SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();
        builder.setPeriodLength(7);
        builder.setSeasonalWidth(5);
        builder.setInnerIterations(1);
        builder.setRobustnessIterations(15);

        SeasonalTrendLoess smoother = builder.buildSmoother(data, fExogenousInputs2);

        SeasonalTrendLoess.Decomposition stl = smoother.decompose();

        double[] seasonal = stl.getSeasonal();
        double[] trendexog = stl.getTrend();
        double[] residual = stl.getResidual();
        double epsilon = 1.0e-2;
        
        assertArrayEquals(fSeas2, seasonal, epsilon);
        assertArrayEquals(fTrendExog2, trendexog, epsilon);
        assertArrayEquals(fResidual2, residual, epsilon);
    }
    
    
	private final double[] fTestData1 = {
	            2.66317301,   7.43055795,  -1.09080058,  -2.0167677,    5.24267195,
	            8.55982888,   6.53102317,   2.21709705,   8.00858286,  13.65079958,
	           10.86605238,   8.40607475,  10.06986866,  20.05654816,  14.45726026,
	           10.89027523,  13.72679558,  22.27027376,  20.969074,    11.84039023,
	           16.04454344,  22.51325541,  22.53214033,  22.76993466,  26.24738251,
	           33.34524443,  26.19829843,  24.79599275,  27.45750402,  35.22640837,
	           29.46536562,  24.90138197,  32.26541659,  37.04771597,  36.61694616,
	           30.39002656,  36.80041998,  41.32473533,  40.51294453,  32.536061,
	           41.3204631,   45.29825622,  40.12113328,  37.02132557,  42.39081772,
	           49.57460472,  45.32171951,  42.62433987,  49.13030534,  53.70515948,
	           49.94818933,  46.57818841,  50.92024186,  59.41632004,  55.68444948,
	           50.40716159,  60.78940733,  63.83491788,  57.77545506,  53.27563911,
	           59.53563549,  64.9965422,   64.25757031,  56.60437994,  63.83775563,
	           68.94140784,  68.09236571,  59.16288794,  67.27500163,  73.75618862,
	           90.63871284,  86.9218058,   91.56842022,  99.97814491,  94.62950756,
	           94.93530211,  92.98335703, 103.24120133,  95.90973492,  92.40398236,
	          103.97016918, 109.4896283,   98.28762904,  97.5544526,  103.8683043,
	          105.73657578, 105.90233898, 102.78668243, 108.43453029, 110.01121246,
	          112.21541647, 106.48908795, 111.87617594, 116.49221408, 115.42391803,
	          111.8365383,  115.03581372, 122.17917523, 119.65399725, 110.09097576  
	};
	
    
    private final double[][] fExogenousInputs1 = {
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
         1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}
    };
	
	private final double[] fSeas1 = {
	         2.62, 6.34, -2.97, -5.03, 1.3, 5.62, -1.35, -4.58, -0.12, 4.52, 0.6, -4.42, -1.29, 5.55, 
	         0.18, -4.24, -1.94, 6.4, 0.69, -5.33, -2.48, 6.48, 0.96, -4.55, -2.12, 6.57, -1.02, -3.43, 
	         -1.48, 5.29, -0.09, -4.29, -0.49, 4.32, 0.96, -5.28, 0.22, 3.68, 1.93, -6.2, 0.99, 4.09, 
	         0.63, -5.81, 1.12, 4.51, -0.73, -5.04, 0.97, 4.91, -0.04, -4.78, -1.42, 5.71, 0.96, -5.41, 
	         -1.07, 5.27, 1.61, -5.98, -0.63, 4.59, 2.33, -6.29, -0.37, 4.23, 1.69, -5.43, -0.67, 4.93, 
	         0.79, -4.64, -1.54, 6.1, 0.1, -4.61, -2.28, 7.4, -0.47, -4.55
	};

    private final double[] fTrendExog1 = {
             0.05, 1.02, 1.99, 2.97, 3.94, 4.96, 5.99, 7.01, 8.04, 9.05, 10.0, 11.03, 12.05, 12.99, 
             13.89, 14.78, 15.51, 15.99, 16.37, 19.42, 21.67, 23.03, 27.4, 26.98, 26.89, 27.13, 27.6, 
             28.22, 29.02, 29.9, 30.78, 31.7, 32.66, 33.61, 34.61, 35.61, 36.58, 37.51, 38.42, 39.32, 
             40.21, 41.11, 42.08, 43.08, 44.08, 45.08, 46.11, 47.13, 48.14, 49.18, 50.22, 51.29, 52.42, 
             53.58, 54.75, 55.82, 56.7, 57.59, 58.48, 59.37, 60.16, 61.02, 61.97, 62.96, 64.0, 65.02, 
             66.01, 66.99, 67.98, 69.05, 90.21, 91.42, 92.55, 93.4, 94.5, 94.81, 95.35, 95.87, 96.39, 96.91
    };
    
    private final double[] fResidual1 = {
            -0.01, 0.07, -0.11, 0.04, -0.0, -2.02, 1.89, -0.21, 0.09, 0.08, 0.26, 1.79, -0.69, 1.51, 0.38, 0.35, 
             0.17, -0.11, 3.9, -2.25, -3.15, -7.0, -5.83, 0.33, 1.47, -0.36, -0.38, 0.01, -0.08, 0.03, -1.22, -2.51, 
             0.1, -0.89, 1.05, 0.06, 0.0, 0.13, 0.16, -0.58, 0.12, 0.1, -2.58, -0.24, -2.81, -0.01, -0.06, 0.54, 0.02, 
            -0.38, -0.23, 0.06, -0.08, 0.13, -0.03, -0.01, 5.15, 0.98, -2.31, -0.12, 0.0, -0.62, -0.04, -0.07, 0.21, 
            -0.31, 0.39, -2.4, -0.03, -0.22, -0.37, 0.15, 0.56, 0.47, 0.03, 4.73, -0.08, -0.03, -0.01, 0.04
    };


    private final double[] fTestData2 = {
            2.66317301,  5.33971536,  1.78383898,  2.4812481 , -2.06058508,
           15.16292013, -3.26524635, -0.71173567,  5.72171146,  3.12345013,
            0.46430326,  1.79448075, -7.93387332, -3.76209373,  1.49148685,
            4.38083085, -0.01178305,  3.55555279, -0.23406307, -9.78712929,
           -7.54373519, 15.28628171,  0.6101824 ,  9.47371641,  3.88544545,
           -2.43696806, -5.84183345, -1.76256471, -4.69758699,  4.04678648,
            4.91126063, -1.4129282 , -1.12729924, -5.83708706, -3.23383794,
            1.14112849,  3.86829227,  1.34750249,  4.45462199, -1.90990723,
           -4.43979393, -5.87917963, -0.39829689,  1.99060575,  1.61409693,
            2.84812744, -3.50237931, -3.25931336, -1.2712182 , -0.2336667 ,
            5.03185603,  8.05394815,  0.52967856, -0.31027063, -1.28135635,
           -2.82539135,  4.87694747,  5.89550188,  5.59732649,  2.38602137,
           -3.8047804 , -4.59135793, -0.24716425, -2.29829135,  5.57201633,
            6.09568361,  3.81265396, -2.59765615, -4.41101586, -5.70579736,
            4.06140176,  5.22263493,  3.13299169,  6.0029356 , -2.23865408,
           -0.72322911, -5.54740987,  0.69159686,  1.2364908 ,  6.34234537,
            7.06776242,  2.3294328 , -8.01877755, -2.90247372, -0.97527576,
           -0.58146171,  7.6236892 ,  3.74787401, -0.94695949, -5.86732334,
           19.26843184,  1.58754629,  7.29028785,  2.68564362,  4.78888686,
            1.73520946, -6.73385209, -2.95223744,  3.37994379, -1.27811177
    };
       
    private final double[][] fExogenousInputs2 = {
       {  0.,  0.,  0.,  0.,  0., 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 1.,  0.,  0.,  0.,  0.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 1.,
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.},
       {       0.        ,  0.        ,  0.        ,  0.3285971 , -1.13383833,
               1.47773081,  0.11288789, -0.92883272,  1.80397118, -1.40198901,
              -2.57116781,  1.5578247 , -1.12910242, -1.90948448,  1.03422659,
              -0.41860179, -2.61321819,  1.11586033, -1.03371838, -2.75287995,
               0.32087879, -1.2269737 , -3.83111534, -0.17085781, -0.53135575,
              -3.61279379, -1.16549232, -0.64940006, -4.15509101, -1.08877931,
               0.57125544, -2.48372887,  0.77670287, -0.01016347, -1.94162669,
               0.75110193, -0.84128511, -2.8518724 , -0.22774123,  1.72345047,
              -0.88561748, -1.26827844,  1.48056983, -0.93987723, -1.65136034,
               1.10410402, -0.65468012,  0.99098633,  1.50763388,  0.06117381,
               1.17450929,  2.60112018, -0.55998199,  0.44282803,  1.90883374,
               0.67660446,  0.08754014,  0.15142659,  0.94723188,  0.94096357,
              -1.17099719,  1.28673943,  1.40442286, -0.90267129,  1.82510328,
               2.27963621, -0.44913045,  2.4088746 ,  1.18862207, -1.55282856,
               3.42268892,  0.39167172, -1.31006808,  1.85537199, -0.69874294,
              -0.78389166,  1.37839051, -0.54960446, -0.58240154,  3.06372344,
               0.92817454,  1.0092232 ,  0.56823297,  1.45223109, -0.84358006,
              -0.2271949 ,  2.84671066,  0.79177288,  0.78792891,  2.99610376,
               0.96217278,  1.09845834,  3.50495449, -0.68121002,  1.19555013,
               2.06808985, -0.89502625,  0.77774475,  1.72594654, -1.27824494}
    };
    
    private final double[] fSeas2 = {
            2.11, 4.7, 1.37, 1.82, -1.37, -4.61, -3.51, 1.23, 4.44, 1.63, 2.18, -0.28, -5.67, 
            -2.93, 0.22, 4.26, 1.94, 2.62, 0.62, -7.37, -2.33, 0.05, 4.43, 2.76, 2.08, -0.9, 
            -6.18, -2.1, 0.06, 4.59, 3.59, 1.1, -2.4, -5.34, -1.45, -0.17, 4.31, 3.78, 1.63, 
            -2.55, -4.85, -2.23, -0.25, 3.79, 3.95, 2.12, -2.48, -4.11, -3.05, -0.09, 3.96, 
            4.36, 1.83, -2.63, -4.3, -3.51, -0.73, 4.48, 4.74, 2.78, -2.81, -4.75, -3.89, 
            -1.46, 4.77, 4.61, 3.68, -2.33, -5.4, -4.48, 0.07, 3.17, 4.02, 4.78, -0.27, -7.06, 
            -4.29, 0.33, 1.62, 3.54, 6.0, 1.26, -8.34, -4.03, 0.55, 0.06
    };

    private final double[] fTrendExog2 = {
            0.48, 0.45, 0.42, 0.64, -0.49, 19.77, 0.43, -0.39, 1.7, -0.76, -1.68, 1.58, -0.53, -1.17, 
            1.22, -0.03, -1.94, 0.97, -0.91, -2.44, 0.45, 15.24, -3.64, 0.69, 0.48, -3.7, -0.28, 0.31, 
            -4.51, -0.57, 1.38, -2.48, 1.33, 0.29, -1.77, 0.72, -0.83, -2.59, -0.54, 0.93, -1.18, -1.55, 
            0.73, -1.51, -2.33, 0.69, -1.16, 0.91, 1.7, -0.22, 1.39, 3.58, -1.25, 0.31, 2.62, 0.75, -0.02, 
            0.13, 0.85, 0.72, -0.9, 0.72, 0.73, -0.69, 1.05, 1.46, -0.37, 1.79, 0.93, -1.12, 2.77, 0.44, 
            -0.92, 1.71, -0.37, -0.41, 1.3, -0.22, -0.34, 2.83, 0.85, 0.85, 0.36, 1.08, -1.14, -0.67
    };
   
    private final double[] fResidual2 = {
            0.07, 0.2, -0.01, 0.02, -0.2, -0.0, -0.18, -1.56, -0.42, 2.26, -0.03, 0.49, -1.74, 0.33, 0.05, 
            0.15, -0.02, -0.04, 0.06, 0.02, -5.66, 0.0, -0.19, 6.03, 1.32, 2.16, 0.62, 0.02, -0.25, 0.03, 
            -0.06, -0.03, -0.05, -0.79, -0.01, 0.59, 0.39, 0.16, 3.37, -0.28, 1.59, -2.1, -0.87, -0.3, -0.0, 
            0.04, 0.13, -0.07, 0.08, 0.08, -0.32, 0.11, -0.05, 2.01, 0.4, -0.07, 5.63, 1.29, 0.0, -1.11, -0.1, 
            -0.56, 2.91, -0.15, -0.25, 0.03, 0.5, -2.06, 0.06, -0.1, 1.22, 1.61, 0.03, -0.48, -1.6, 6.75, -2.56, 
            0.58, -0.04, -0.03, 0.22, 0.22, -0.03, 0.05, -0.38, 0.03
    };


}
