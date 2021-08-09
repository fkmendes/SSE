package test;

import SSE.ConstantLinkFn;
import SSE.NormalCenteredAtObservedLinkFn;
import SSE.QuaSSEDistribution;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import SSE.LogisticFunction;
import beast.core.parameter.RealParameter;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;

public class QuaSSEFunctionsTest {

    final static Double EPSILON = 1e-6;

    static Double[] x0, y1, y0, r;
    static double dx;
    static RealParameter dtrp, diffusionrp;
    double[] xRuler, lambdaHi, muHi;

    static LogisticFunction lfn;
    static ConstantLinkFn cfn;

    /*
     * Applies logistic function to many x values, checks return
     *
     * We're only looking at high resolution here.
     */
    @Test
    public void testLogistic() {

        int[] nUsefulBins = new int[] { 999, 3999 }; // lo, hi
        double x2Add; // x0 is xmid, x2Add is xmin
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

        lfn = new LogisticFunction();
        lfn.initByName( "curveMaxBase", y0rp, "added2CurveMax", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);
        double[] lfnOut = lfn.getY(xRuler, lambdaHi, true);

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

        cfn = new ConstantLinkFn();
        cfn.initByName("yV", yValuerp);
        double[] cfnOut = cfn.getY(xRuler, lambdaHi, true);

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
        double[] lfnOut = lfn.getY(xRuler, lambdaHi, true);
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

        ConstantLinkFn cfn = new ConstantLinkFn();
        cfn.initByName("yV", yValuerp);
        double[] cfnOut = cfn.getY(xRuler, muHi, true);
    }
}
