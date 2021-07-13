package test;

import SSE.ConstantFunction;
import SSE.QuaSSEDistribution;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import SSE.LogisticFunction;
import beast.core.parameter.RealParameter;
import org.junit.Assert;
import org.junit.Test;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

public class QuaSSEDimensionsAndFunctionsTest {

    final static Double EPSILON = 1e-6;

    QuaSSEDistribution q;

    Double[] x0, y1, y0, r;
    double dx;
    double[] xRuler, lambdaHi, muHi;

    LogisticFunction lfn;
    ConstantFunction cfn;

    /*
     * Applies logistic function to many x values, checks return
     *
     * We're only looking at high resolution here.
     */
    @Test
    public void testLogistic() {

        int[] nUsefulBins = new int[] { 999, 3999 }; // lo, hi
        double dx, x2Add; // x0 is xmid, x2Add is xmin
        xRuler = new double[nUsefulBins[1]];
        lambdaHi = new double[nUsefulBins[1]];
        x2Add = -4.9975;
        dx = 0.01;
        for (int i=0; i<nUsefulBins[1]; i++) {
            xRuler[i] = x2Add;
            x2Add += (dx / 4.0); // high-res : low-res
        }

        // RealParameter x = new RealParameter(xRuler);

        // logistic realparameter's
        x0 = new Double[] { 0.0 };
        y1 = new Double[] { 0.2 };
        y0 = new Double[] { 0.1 };
        r = new Double[] { 2.5 }; // logistic function parameters
        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);

        LogisticFunction lfn = new LogisticFunction();
        lfn.initByName( "curveMaxBase", y0rp, "added2CurveMax", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);
        double[] lfnOut = lfn.getMacroParams(xRuler, lambdaHi);

        double[] expectedLfnOut1to10 = new double[] { 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004 };
        double[] expectedLfnOut2500to2510 = new double[] { 0.195816352937447, 0.195841335181637, 0.195866174682834, 0.195890872179954, 0.195915428409005, 0.195939844103087, 0.19596411999239, 0.195988256804195, 0.196012255262877, 0.196036116089903 };
        double[] expectedLfnOut3989to3999 = new double[] { 0.199999600814288, 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637};

        assertArrayEquals(expectedLfnOut1to10, Arrays.copyOfRange(lfnOut, 0, 10), EPSILON);
        assertArrayEquals(expectedLfnOut2500to2510, Arrays.copyOfRange(lfnOut, 2500, 2510), EPSILON);
        assertArrayEquals(expectedLfnOut3989to3999, Arrays.copyOfRange(lfnOut, 3988, 3999), EPSILON);
    }

    /*
     * Applies constant function to many x values, checks return
     *
     * We're only looking at high resolution here.
     */
    @Test
    public void testConstant() {
        int[] nUsefulBins = new int[] { 999, 3999 }; // lo, hi
        xRuler = new double[nUsefulBins[1]];
        lambdaHi = new double[nUsefulBins[1]];

        double x2Add = -4.9975;
        dx = 0.01;
        for (int i=0; i<nUsefulBins[1]; i++) {
            xRuler[i] = x2Add;
            x2Add += (dx / 4.0); // high-res : low-res
        }

        // constant realparameter's
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValuerp = new RealParameter(yValue);

        cfn = new ConstantFunction();
        cfn.initByName("yV", yValuerp);
        double[] cfnOut = cfn.getMacroParams(xRuler, lambdaHi);

        double[] expectedCfn = new double[3999];
        for (int i=0; i<expectedCfn.length; i++) {
            expectedCfn[i] = 0.03;
        }

        assertArrayEquals(expectedCfn, cfnOut, 0.0);
    }

    /*
     * If we try to get the macroevolutionary parameters after applying
     * the linking function (e.g., logistic) -- but forgetting to set
     * the quantitative trait array (ruler) and the size of y (i.e.,
     * the macroevol parameter array) -- we get an error.
     *
     * This setting is done by the QuaSSE likelihood in initialization.
     */
    @Test(expected = RuntimeException.class)
    public void testQu2MacroevolFailLogistic() {

        int nUsefulBins = 3999;
        int someWrongNumber = 3998;
        double dx, x2Add; // x0 is xmid, x2Add is xmin
        xRuler = new double[nUsefulBins];
        lambdaHi = new double[someWrongNumber];
        x2Add = -4.9975;
        dx = 0.01;
        for (int i=0; i<nUsefulBins; i++) {
            xRuler[i] = x2Add;
            x2Add += (dx / 4.0); // high-res : low-res
        }

        // logistic realparameter's
        Double[] x0, y1, y0, r;
        x0 = new Double[] { 0.0 };
        y1 = new Double[] { 0.2 };
        y0 = new Double[] { 0.1 };
        r = new Double[] { 2.5 }; // logistic function parameters
        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);

        LogisticFunction lfn = new LogisticFunction();
        lfn.initByName( "curveMaxBase", y0rp, "added2CurveMax", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);
        double[] lfnOut = lfn.getMacroParams(xRuler, lambdaHi);
    }

    /*
     * If we try to get the macroevolutionary parameters after applying
     * the linking function (e.g., constant) -- but forgetting to set
     * the quantitative trait array (ruler) and the size of y (i.e.,
     * the macroevol parameter array) -- we get an error.
     *
     * This setting is done by the QuaSSE likelihood in initialization.
     */
    @Test(expected = RuntimeException.class)
    public void testQu2MacroevolFailConstant() {

        int nUsefulBins = 3999;
        int someWrongNumber = 3998;
        double dx, x2Add; // x0 is xmid, x2Add is xmin
        xRuler = new double[nUsefulBins];
        muHi = new double[someWrongNumber];
        x2Add = -4.9975;
        dx = 0.01;
        for (int i=0; i<nUsefulBins; i++) {
            xRuler[i] = x2Add;
            x2Add += (dx / 4.0); // high-res : low-res
        }

        // constant realparameter's
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValuerp = new RealParameter(yValue);

        ConstantFunction cfn = new ConstantFunction();
        cfn.initByName("yV", yValuerp);
        double[] cfnOut = cfn.getMacroParams(xRuler, muHi);
    }

    /*
     * Checks that QuaSSE likelihood class it setting its own
     * dimensions correctly, and then applying a logistic function
     * on the quantitative trait for lambda, and a constant function
     * for mu.
     */
    @Test
    public void testDimensions() {

        // tree
        String treeStr = "((sp1:1.0,sp2:1.0):1.0,sp3:2.0);";
        Tree myTree = new TreeParser(treeStr, false, false, true, 0);

        // logistic realparameter's for lambda
        x0 = new Double[] { 0.0 };
        y1 = new Double[] { 0.2 };
        y0 = new Double[] { 0.1 };
        r = new Double[] { 2.5 }; // logistic function parameters
        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);

        lfn = new LogisticFunction();
        lfn.initByName( "curveMaxBase", y0rp, "added2CurveMax", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);

        // constant realparameter's for mu
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValuerp = new RealParameter(yValue);

        ConstantFunction cfn = new ConstantFunction();
        cfn.initByName("yV", yValuerp);

        // dimension stuff
        Double[] dt = new Double[] { 0.05 };
        RealParameter dtrp = new RealParameter(dt);

        Double[] xMid = new Double[] { 0.0 };
        RealParameter xMidrp = new RealParameter(xMid);

        Double[] dxBin = new Double[] { 0.01 };
        RealParameter dxBinrp = new RealParameter(dxBin);

        Double[] flankWidthScaler = new Double[] { 5.0 };
        RealParameter flankWidthScalerrp = new RealParameter(flankWidthScaler);

        Integer[] hiLoRatio = new Integer[] { 4 };
        IntegerParameter hiLoRatiorp = new IntegerParameter(hiLoRatio);

        Integer[] nXbins = new Integer[] { 1024 };
        IntegerParameter nXbinsip = new IntegerParameter(nXbins);

        // qu trait stuff
        Double[] drift = new Double[] { 0.0 };
        RealParameter driftrp = new RealParameter(drift);

        Double[] diffusion = new Double[] { 0.01 };
        RealParameter diffusionrp = new RealParameter(diffusion);

        // QuaSSE stuff
        q = new QuaSSEDistribution();
        q.initByName("dt", dtrp, "nX", nXbinsip, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftrp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree);

        int nXbinsLo = q.getnXbins(true);
        int nXbinsHi = q.getnXbins(false);
        int nUsefulXbinsLo = q.getNUsefulXbins(true);
        int nUsefulXbinsHi = q.getNUsefulXbins(false);
        int nLeftFlanksLo = q.getNLeftFlanks(true);
        int nLeftFlanksHi = q.getNLeftFlanks(false);
        int nRightFlanksLo = q.getNRightFlanks(true);
        int nRightFlanksHi = q.getNRightFlanks(false);
        double xMinLo = q.getXMinLo();
        double xMinHi = q.getXMinHi();
        double[] xLo = q.getX(true);
        double[] xHi = q.getX(false);
        double[] lambdaLo = q.getLambda(true);
        double[] lambdaHi = q.getLambda(false);
        double[] muLo = q.getMu(true);
        double[] muHi = q.getMu(false);

        // expected x rulers
        double[] expectedXLoFirst10 = new double[] { -4.99, -4.98, -4.97, -4.96, -4.95, -4.94, -4.93, -4.92, -4.91, -4.90 };
        double[] expectedXLoLast10 = new double[] { 4.90, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99 };
        double[] expectedXHiFirst10 = new double[] { -4.9975, -4.995, -4.9925, -4.99, -4.9875, -4.985, -4.9825, -4.98, -4.9775, -4.975 };
        double[] expectedXHiLast10 = new double[] { 4.975, 4.9775, 4.98, 4.9825, 4.985, 4.9875, 4.99, 4.9925, 4.995, 4.9975 };

        // expected lambdas
        double[] expectedLambdaLoFirt10 = new double[] { 0.100000382097925, 0.100000391770742, 0.100000401688425, 0.100000411857174, 0.100000422283344, 0.100000432973452, 0.100000443934178, 0.100000455172374, 0.100000466695064, 0.10000047850945 };
        double[] expectedLambdaLoLast10 = new double[] { 0.1999995, 0.1999995, 0.1999995, 0.1999996, 0.1999996, 0.1999996, 0.1999996, 0.1999996, 0.1999996, 0.1999996 };
        double[] expectedLambdaHiFirst10 = new double[] { 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004, 0.1000004 };
        double[] expectedLambdaHiLast10 = new double[] { 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637 };

        // expected mus
        double[] expectedMuLoHiFirstLast10 = new double[] { 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03 };

        Assert.assertEquals(nXbinsLo, 1024, 0.0);
        Assert.assertEquals(nXbinsHi, 1024 * 4, 0.0);
        Assert.assertEquals(999, nUsefulXbinsLo, 0.0);
        Assert.assertEquals(3999, nUsefulXbinsHi, 0.0);
        Assert.assertEquals(12, nLeftFlanksLo, 0.0);
        Assert.assertEquals(48, nLeftFlanksHi, 0.0);
        Assert.assertEquals(12, nRightFlanksLo, 0.0);
        Assert.assertEquals(48, nRightFlanksHi, 0.0);
        Assert.assertEquals(-4.99, xMinLo, EPSILON);
        Assert.assertEquals(-4.9975, xMinHi, EPSILON);
        Assert.assertEquals(999, xLo.length, 0.0);
        Assert.assertEquals(3999, xHi.length, 0.0);
        Assert.assertArrayEquals(expectedXLoFirst10, Arrays.copyOfRange(xLo, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedXLoLast10, Arrays.copyOfRange(xLo, 989, 999), EPSILON);
        Assert.assertArrayEquals(expectedXHiFirst10, Arrays.copyOfRange(xHi, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedXHiLast10, Arrays.copyOfRange(xHi, 3989, 3999), EPSILON);
        Assert.assertEquals(999, lambdaLo.length, 0.0);
        Assert.assertEquals(3999, lambdaHi.length, 0.0);
        Assert.assertArrayEquals(expectedLambdaLoFirt10, Arrays.copyOfRange(lambdaLo, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedLambdaLoLast10, Arrays.copyOfRange(lambdaLo, 989, 999), EPSILON);
        Assert.assertArrayEquals(expectedLambdaHiFirst10, Arrays.copyOfRange(lambdaHi, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedLambdaHiLast10, Arrays.copyOfRange(lambdaHi, 3989, 3999), EPSILON);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 0, 10), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 989, 999), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 0, 10), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 3989, 3999), 0.0);
    }

}
