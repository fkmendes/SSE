package test;

import SSE.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static SSE.SSEUtils.everyOtherInPlace;
import static SSE.SSEUtils.propagateEandDinTQuaSSEInPlace;

public class QuaSSEDistributionTest {

    final static Double EPSILON = 1e-6;
    final static Double EPSILON2 = 1e-16;
    final static Double EPSILON3 = 1e-14;

    static QuaSSEDistribution q, q2, q3;
    static Tree myTree, myTree2;
    static int nDimensionsE, nDimensionsD;
    double[] birthRate, deathRate;

    static List<Double> data, data2;
    static RealParameter quTraitrp, quTraitrp2;
    static RealParameter dtrp, tcrp, diffusionrp;

    static Double[] x0, y1, y0, r;
    static LogisticFunction lfn;
    static ConstantLinkFn cfn;
    static NormalCenteredAtObservedLinkFn nfn, nfn2;

    @BeforeClass
    public static void setupQuaSSEDist() {
        // tree
        String treeStr = "((sp1:0.5787065,sp2:0.5787065):1.0002893,sp3:1.5789958);";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        String treeStr2 = "(sp1:0.01,sp2:0.01);";
        myTree2 = new TreeParser(treeStr2, false, false, true, 0);

        // qu trait data
        String spNames = "sp1 sp2 sp3";
        data = Arrays.asList(-0.061675531, 0.001079479, 0.082005868);
        quTraitrp = new RealParameter();
        quTraitrp.initByName("value", data, "keys", spNames);

        String spNames2 = "sp1 sp2";
        data2 = Arrays.asList(0.0, 0.1);
        quTraitrp2 = new RealParameter();
        quTraitrp2.initByName("value", data2, "keys", spNames2);

        // qu trait stuff
        Double[] drift = new Double[] { 0.0 };
        RealParameter driftrp = new RealParameter(drift);

        Double[] diffusion = new Double[] { 0.001 };
        diffusionrp = new RealParameter(diffusion);

        // dimension stuff
        Double[] dt = new Double[] { 0.01 };
        dtrp = new RealParameter(dt);

        Double[] tc = new Double[] { 100.0 };
        tcrp = new RealParameter(tc);

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

        cfn = new ConstantLinkFn();
        cfn.initByName("yV", yValuerp);

        // link function for D's

        Double[] sdNormaQuTraitValue = new Double[] { 0.05 };
        RealParameter sdNormaQuTraitValuerp = new RealParameter(sdNormaQuTraitValue);
        nfn = new NormalCenteredAtObservedLinkFn();
        nfn.initByName("quTraits", quTraitrp, "sdNormalQuTrValue", sdNormaQuTraitValuerp);
        nfn2 = new NormalCenteredAtObservedLinkFn();
        nfn2.initByName("quTraits", quTraitrp2, "sdNormalQuTrValue", sdNormaQuTraitValuerp);

        Double[] xMid = new Double[] { 0.0 };
        RealParameter xMidrp = new RealParameter(xMid);

        Double[] dxBin = new Double[] { 0.01 };
        RealParameter dxBinrp = new RealParameter(dxBin);

        Double[] flankWidthScaler = new Double[] { 10.0 };
        RealParameter flankWidthScalerrp = new RealParameter(flankWidthScaler);

        Integer[] hiLoRatio = new Integer[] { 4 };
        IntegerParameter hiLoRatiorp = new IntegerParameter(hiLoRatio);

        Integer[] nXbins = new Integer[] { 1024 };
        IntegerParameter nXbinsip = new IntegerParameter(nXbins);

        Integer[] nXbins2 = new Integer[] { 48 };
        IntegerParameter nXbinsip2 = new IntegerParameter(nXbins2);

        // QuaSSE stuff

        // more bins and 3-sp tree (more complex)
        q = new QuaSSEDistribution();
        q.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbinsip, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftrp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree,
                "q2d", nfn);

        // fewer bins and 2-sp tree (less complex)
        q2 = new QuaSSEDistribution();
        q2.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbinsip2, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftrp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree2,
                "q2d", nfn2);

        // same as q2, but needs a new instance, otherwise its starting values will have been processed by previous unit tests
        q3 = new QuaSSEDistribution();
        q3.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbinsip2, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftrp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree2,
                "q2d", nfn2);
    }

    /*
     * Checks that QuaSSE likelihood class it setting its own
     * dimensions correctly, and then applying a logistic function
     * on the quantitative trait for lambda, and a constant function
     * for mu.
     */
//    @Test
//    public void testDimensions() {
//
//        int nXbinsLo = q.getnXbins(true);
//        int nXbinsHi = q.getnXbins(false);
//        int nUsefulXbinsLo = q.getNUsefulXbins(true);
//        int nUsefulXbinsHi = q.getNUsefulXbins(false);
//        int nLeftFlanksLo = q.getNLeftFlanks(true);
//        int nLeftFlanksHi = q.getNLeftFlanks(false);
//        int nRightFlanksLo = q.getNRightFlanks(true);
//        int nRightFlanksHi = q.getNRightFlanks(false);
//        double xMinLo = q.getXMinLo();
//        double xMinHi = q.getXMinHi();
//        double[] xLo = q.getX(true);
//        double[] xHi = q.getX(false);
//        double[] lambdaLo = q.getLambda(true);
//        double[] lambdaHi = q.getLambda(false);
//        double[] muLo = q.getMu(true);
//        double[] muHi = q.getMu(false);
//
//        // expected x rulers
//        double[] expectedXLoFirst10 = new double[] { -4.99, -4.98, -4.97, -4.96, -4.95, -4.94, -4.93, -4.92, -4.91, -4.90 };
//        double[] expectedXLoLast10 = new double[] { 4.90, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99 };
//        double[] expectedXHiFirst10 = new double[] { -4.9975, -4.995, -4.9925, -4.99, -4.9875, -4.985, -4.9825, -4.98, -4.9775, -4.975 };
//        double[] expectedXHiLast10 = new double[] { 4.975, 4.9775, 4.98, 4.9825, 4.985, 4.9875, 4.99, 4.9925, 4.995, 4.9975 };
//
//        // expected lambdas
//        double[] expectedLambdaLoFirt10 = new double[] { 0.100000382097925, 0.100000391770742, 0.100000401688425, 0.100000411857174, 0.100000422283344, 0.100000432973452, 0.100000443934178, 0.100000455172374, 0.100000466695064, 0.10000047850945 };
//        double[] expectedLambdaLoLast10 = new double[] { 0.199999521490551, 0.199999533304936, 0.199999544827626, 0.199999556065822, 0.199999567026548, 0.199999577716656, 0.199999588142826, 0.199999598311575, 0.199999608229258, 0.199999617902075 };
//        double[] expectedLambdaHiFirst10 = new double[] { 0.100000375000363, 0.100000377351446, 0.100000379717269, 0.100000382097925, 0.100000384493506, 0.100000386904106, 0.10000038932982, 0.100000391770742, 0.100000394226967, 0.100000396698591 };
//        double[] expectedLambdaHiLast10 = new double[] { 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637 };
//
//        // expected mus
//        double[] expectedMuLoHiFirstLast10 = new double[] { 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03 };
//
//        Assert.assertEquals(nXbinsLo, 1024, 0.0);
//        Assert.assertEquals(nXbinsHi, 1024 * 4, 0.0);
//        Assert.assertEquals(999, nUsefulXbinsLo, 0.0);
//        Assert.assertEquals(3999, nUsefulXbinsHi, 0.0);
//        Assert.assertEquals(12, nLeftFlanksLo, 0.0);
//        Assert.assertEquals(48, nLeftFlanksHi, 0.0);
//        Assert.assertEquals(12, nRightFlanksLo, 0.0);
//        Assert.assertEquals(48, nRightFlanksHi, 0.0);
//        Assert.assertEquals(-4.99, xMinLo, EPSILON);
//        Assert.assertEquals(-4.9975, xMinHi, EPSILON);
//        Assert.assertEquals(999, xLo.length, 0.0);
//        Assert.assertEquals(3999, xHi.length, 0.0);
//        Assert.assertArrayEquals(expectedXLoFirst10, Arrays.copyOfRange(xLo, 0, 10), EPSILON);
//        Assert.assertArrayEquals(expectedXLoLast10, Arrays.copyOfRange(xLo, 989, 999), EPSILON);
//        Assert.assertArrayEquals(expectedXHiFirst10, Arrays.copyOfRange(xHi, 0, 10), EPSILON);
//        Assert.assertArrayEquals(expectedXHiLast10, Arrays.copyOfRange(xHi, 3989, 3999), EPSILON);
//        Assert.assertEquals(999, lambdaLo.length, 0.0);
//        Assert.assertEquals(3999, lambdaHi.length, 0.0);
//        Assert.assertArrayEquals(expectedLambdaLoFirt10, Arrays.copyOfRange(lambdaLo, 0, 10), EPSILON);
//        Assert.assertArrayEquals(expectedLambdaLoLast10, Arrays.copyOfRange(lambdaLo, 989, 999), EPSILON);
//        Assert.assertArrayEquals(expectedLambdaHiFirst10, Arrays.copyOfRange(lambdaHi, 0, 10), EPSILON);
//        Assert.assertArrayEquals(expectedLambdaHiLast10, Arrays.copyOfRange(lambdaHi, 3989, 3999), EPSILON);
//        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 0, 10), 0.0);
//        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 989, 999), 0.0);
//        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 0, 10), 0.0);
//        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 3989, 3999), 0.0);
//    }

    /*
     * Checks that three species in a 3sp-tree get their D's correctly
     * initialized (using standard QuaSSE initialization)
     */
    @Test
    public void testInitializationOfTips() {
        double[][][] esDsHi = q.getEsDs(false);

        double[] expectedSp1Ds = new double[] { 8.13623712450468e-08, 1.1005578408194e-07, 1.48496566106919e-07, 1.99863833658288e-07, 2.68328176840003e-07, 3.59345830399246e-07, 4.80035331999228e-07, 6.39658332527801e-07, 8.50231484597812e-07, 1.12730275515906e-06 };
        double[] expectedSp2Ds = new double[] { 1.82648393304909e-11, 2.63062707587468e-11, 3.77934883411903e-11, 5.41612820961595e-11, 7.74239202284714e-11, 1.10401669802924e-10, 1.57032805818765e-10, 2.2280216312045e-10, 3.15328103714112e-10, 4.45164187528514e-10 };
        double[] expectedSp3Ds = new double[] { 3.51799453441839e-17, 5.49394711042998e-17, 8.55831075496442e-17, 1.32985990766802e-16, 2.06128478871848e-16, 3.18701690685798e-16, 4.91524307682648e-16, 7.56170789433044e-16, 1.16040357229278e-15, 1.77628431652675e-15 };

        Assert.assertArrayEquals(expectedSp1Ds, Arrays.copyOfRange(esDsHi[0][1], 1885, 1895), EPSILON);
        Assert.assertArrayEquals(expectedSp2Ds, Arrays.copyOfRange(esDsHi[1][1], 1885, 1895), EPSILON);
        Assert.assertArrayEquals(expectedSp3Ds, Arrays.copyOfRange(esDsHi[2][1], 1885, 1895), EPSILON);
    }

    /*
     * Checks that propagate methods in time for E's and D's
     * inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 48 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes48BinsOutsideClassJustT() {

        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q2.getEsDsAtNode(nodeIdx, true);
        // printing
        // System.out.println("D's before prop in t: " + Arrays.toString(esDsHiAtNode[1])); // D's

        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i=0; i < esDsLoAtNode[ithDim].length; i++) {
                scratchAtNode[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }

        // just propagate in t, in place
        q2.propagateTInPlace(esDsLoAtNode, scratchAtNode, true);

        esDsLoAtNode = q2.getEsDsAtNode(nodeIdx, true);
        double[] esHiAtNode = esDsLoAtNode[0];
        double[] dsHiAtNode = esDsLoAtNode[1];

        // printing
        //  for (int i=0; i<esDsHiAtNode[0].length; ++i) {
        //      System.out.println("e" + i + " = " + esHiAtNode[i] + " d" + i + " = " + dsHiAtNode[i]);
        //  }

        double[] expectedSp1EsAfterPropT = new double[] { 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropT = new double[] { 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esHiAtNode, 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(dsHiAtNode, 0, 48), EPSILON3);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 48 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes48BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsHiAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsHiAtNode = q3.getEsDsAtNode(nodeIdx, true);

        // making deep copy for assert below (checking that initial D's are correct)
        double[][] esDsHiAtNodeInitial = new double[esDsHiAtNode.length][esDsHiAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsHiAtNode[0], esDsHiAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsHiAtNode[1], esDsHiAtNode[1].length); // D

        // just propagate in x, in place
        double[] fftedfY = q3.getfY(true);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q3.getnXbins(true),0, 0, 1.0); // getting real part for assert below

        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsHiAtNode[0].length];
        q3.propagateXInPlace(esDsHiAtNode, scratchAtNode, true);

        esDsHiAtNode = q3.getEsDsAtNode(nodeIdx, true); // TODO: something wrong, probably with how scratch is being used... need to look at each individual output of fftR.propagate.x and compare to mine
        double[] esLoAtNode = esDsHiAtNode[0];
        double[] dsLoAtNode = esDsHiAtNode[1];

        double[] expectedInitialDs = new double[] { 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113, 0.990144725007753, 0.988438851882212, 0.986703282961364, 0.984967714317709, 0.983261842004857, 0.981614853950298, 0.980054930543287, 0.978608762462816, 0.977301093995625, 0.976154299658014, 0.97518800136532, 0.97441873269917, 0.97385965601673, 0.973520337242051, 0.973406582192705, 0.973520337242051, 0.97385965601673, 0.97441873269917, 0.97518800136532, 0.976154299658014, 0.977301093995625, 0.978608762462816, 0.980054930543287, 0.981614853950298, 0.983261842004857, 0.984967714317709, 0.986703282961364, 0.988438851882212, 0.990144725007753, 0.991791714355113, 0.993351639446935, 0.994797809489411, 0.996105480062091, 0.997252276505192, 0.998218576759891, 0.998987847110852, 0.999546925086093, 0.999886244673461 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropX = new double[] { 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0888278883396178, 0.158599420774613, 0.272077439492024, 0.44845817681081, 0.710214708266689, 1.0806760140649, 1.57993546382425, 2.21932565164129, 2.99530083804265, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661662, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.08067601406491, 0.710214708266689, 0.448458176810809, 0.272077439492024, 0.158599420774613, 0.088827888339617, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(realFFTedfY, 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), EPSILON3);
        Assert.assertArrayEquals(expectedSp1DsAfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), EPSILON3);
    }

    // now we'll test the integration inside the class
    // double[][][] esDsHi = q2.getEsDs(false);

    // TODO: in another test, do inside class
    // q2.processBranch(myTree2.getNode(nodeIdx));
}
