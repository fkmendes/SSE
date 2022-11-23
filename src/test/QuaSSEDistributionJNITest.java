import SSE.ConstantLinkFn;
import SSE.LogisticFunction;
import SSE.NormalCenteredAtObservedLinkFn;
import SSE.QuaSSEDistribution;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

/**
 * @author Kylie Chen
 * @author Fabio K Mendes
 */
public class QuaSSEDistributionJNITest {

    static QuaSSEDistribution q32Dt001;

    static RealParameter dt001Rp;
    static BooleanParameter dynDtbpTrue;
    static RealParameter tc100Rp;
    static IntegerParameter nXbins32Ip;
    static RealParameter dxBin001Rp;
    static RealParameter xMid00Rp;
    static RealParameter flankWidthScaler10Rp;
    static IntegerParameter hiLoRatioIp;
    static RealParameter driftRp;
    static RealParameter diffusionRp0001;
    static LogisticFunction lfn;
    static ConstantLinkFn cfn;
    static NormalCenteredAtObservedLinkFn nfn2Sp;
    static String rootPriorType;

    static Tree bifTreeHeight001;
    static List<Double> data2Sp;
    static RealParameter quTrait2SpRp;

    static Double[] x0, y1, y0, r;

    @BeforeClass
    public static void setupParameters() {
        // tree
        String bifTreeStr001 = "(sp1:0.01,sp2:0.01);";
        bifTreeHeight001 = new TreeParser(bifTreeStr001, false, false, true, 0);

        // qu trait data
        String spNames2Sp = "sp1 sp2";
        data2Sp = Arrays.asList(0.0, 0.1);
        quTrait2SpRp = new RealParameter();
        quTrait2SpRp.initByName("value", data2Sp, "keys", spNames2Sp);

        // dimensions
        Double[] dtD001 = new Double[] { 0.01 };
        dt001Rp = new RealParameter(dtD001);
        // adjust dt dynamically to maximize accuracy and match diversitree
        Boolean[] dynDtTrue = new Boolean[] { true };
        dynDtbpTrue = new BooleanParameter(dynDtTrue);
        Double[] tc100 = new Double[] { 100.0 };
        tc100Rp = new RealParameter(tc100);

        // bins
        Integer[] nXbins32 = new Integer[] { 32 };
        nXbins32Ip = new IntegerParameter(nXbins32);
        Double[] dxBin001 = new Double[] { 0.01 };
        dxBin001Rp = new RealParameter(dxBin001);

        // link functions
        Double[] xMid00 = new Double[] { 0.0 };
        xMid00Rp = new RealParameter(xMid00);
        Double[] flankWidthScaler10 = new Double[] { 10.0 };
        flankWidthScaler10Rp = new RealParameter(flankWidthScaler10);
        Integer[] hiLoRatio = new Integer[] { 4 };
        hiLoRatioIp = new IntegerParameter(hiLoRatio);

        // qu trait stuff
        Double[] drift = new Double[] { 0.0 };
        driftRp = new RealParameter(drift);
        Double[] diffusion0001 = new Double[] { 0.001 };
        diffusionRp0001 = new RealParameter(diffusion0001);

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
        lfn.initByName( "curveYBaseValue", y0rp, "curveMaxY", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);

        // constant realparameter's for mu
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValuerp = new RealParameter(yValue);

        cfn = new ConstantLinkFn();
        cfn.initByName("yV", yValuerp);

        // link funcs for D's
        Double[] sdNormaQuTraitValue005 = new Double[] { 0.05 };
        RealParameter sdNormaQuTraitValue005Rp = new RealParameter(sdNormaQuTraitValue005);
        nfn2Sp = new NormalCenteredAtObservedLinkFn();
        nfn2Sp.initByName("quTraits", quTrait2SpRp, "sdNormalQuTrValue", sdNormaQuTraitValue005Rp);

        // root prior
        rootPriorType = "Observed";
    }

    @Before
    public void setupSSELiks() {
        q32Dt001 = new QuaSSEDistribution();
        q32Dt001.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);
    }

    @Test
    public void testQuaSSEStashJNIInitialization() {
        System.out.println(q32Dt001.getNUsefulTraitBins(true));
        assert(true);
    }


}
