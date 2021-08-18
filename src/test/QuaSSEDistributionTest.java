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

/*
 * Tests for the main methods in QuaSSEDistribution.
 *
 * Note that expectations for D's in
 *
 * (i) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX
 * (ii) testIntegrateOneBranchLoRes4096BinsOutsideClassJustX
 *
 * come from R running on an 2014 iMac with Mac OS X Catalina.
 * These might not be met if tests are being run from machines with
 * different CPU architectures (and small differences are expected
 * even within the same architecture
 */
public class QuaSSEDistributionTest {

    final static Double EPSILON = 1e-6;
    final static Double EPSILON2 = 1e-16;
    final static Double EPSILON3 = 1e-14;

    static QuaSSEDistribution q1024, q48One, q48Two;
    static Tree myTree;
    static int nDimensionsE, nDimensionsD;
    double[] birthRate, deathRate;

    static List<Double> data;
    static RealParameter quTraitrp, quTraitRp;
    static RealParameter dtrp, tcrp, diffusionrp;

    static Double[] x0, y1, y0, r;
    static LogisticFunction lfn;
    static ConstantLinkFn cfn;
    static NormalCenteredAtObservedLinkFn nfn;

    @BeforeClass
    public static void setupQuaSSEDist() {
        // tree
        String treeStr = "(sp1:0.01,sp2:0.01);";
        myTree = new TreeParser(treeStr, false, false, true, 0);

        // qu trait data
        String spNames = "sp1 sp2";
        data = Arrays.asList(0.0, 0.1);
        quTraitRp = new RealParameter();
        quTraitRp.initByName("value", data, "keys", spNames);

        // qu trait stuff
        Double[] drift = new Double[] { 0.0 };
        RealParameter driftRp = new RealParameter(drift);

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
        nfn.initByName("quTraits", quTraitRp, "sdNormalQuTrValue", sdNormaQuTraitValuerp);
        // nfn2 = new NormalCenteredAtObservedLinkFn();
        // nfn2.initByName("quTraits", quTraitRp, "sdNormalQuTrValue", sdNormaQuTraitValuerp);

        Double[] xMid = new Double[] { 0.0 };
        RealParameter xMidrp = new RealParameter(xMid);

        Double[] dxBin = new Double[] { 0.01 };
        RealParameter dxBinrp = new RealParameter(dxBin);

        Double[] flankWidthScaler = new Double[] { 10.0 };
        RealParameter flankWidthScalerrp = new RealParameter(flankWidthScaler);

        Integer[] hiLoRatio = new Integer[] { 4 };
        IntegerParameter hiLoRatiorp = new IntegerParameter(hiLoRatio);

        Integer[] nXbins1024 = new Integer[] { 1024 };
        IntegerParameter nXbins1024Ip = new IntegerParameter(nXbins1024);

        Integer[] nXbins2 = new Integer[] { 48 };
        IntegerParameter nXbins48Ip = new IntegerParameter(nXbins2);

        // QuaSSE stuff

        // more bins and 3-sp tree (more complex)
        q1024 = new QuaSSEDistribution();
        q1024.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbins1024Ip, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftRp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree,
                "q2d", nfn);

        // fewer bins and 2-sp tree (less complex)
        q48One = new QuaSSEDistribution();
        q48One.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbins48Ip, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftRp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree,
                "q2d", nfn);

        // same as q2, but needs a new instance, otherwise its starting values will have been processed by previous unit tests
        q48Two = new QuaSSEDistribution();
        q48Two.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbins48Ip, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftRp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree,
                "q2d", nfn);
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
        double[][][] esDsHi = q1024.getEsDs(false);

        // TODO: update it to simpler tree with just two species
//        double[] expectedSp1Ds = new double[] { 8.13623712450468e-08, 1.1005578408194e-07, 1.48496566106919e-07, 1.99863833658288e-07, 2.68328176840003e-07, 3.59345830399246e-07, 4.80035331999228e-07, 6.39658332527801e-07, 8.50231484597812e-07, 1.12730275515906e-06 };
//        double[] expectedSp2Ds = new double[] { 1.82648393304909e-11, 2.63062707587468e-11, 3.77934883411903e-11, 5.41612820961595e-11, 7.74239202284714e-11, 1.10401669802924e-10, 1.57032805818765e-10, 2.2280216312045e-10, 3.15328103714112e-10, 4.45164187528514e-10 };
//        double[] expectedSp3Ds = new double[] { 3.51799453441839e-17, 5.49394711042998e-17, 8.55831075496442e-17, 1.32985990766802e-16, 2.06128478871848e-16, 3.18701690685798e-16, 4.91524307682648e-16, 7.56170789433044e-16, 1.16040357229278e-15, 1.77628431652675e-15 };
//
//        Assert.assertArrayEquals(expectedSp1Ds, Arrays.copyOfRange(esDsHi[0][1], 1885, 1895), EPSILON);
//        Assert.assertArrayEquals(expectedSp2Ds, Arrays.copyOfRange(esDsHi[1][1], 1885, 1895), EPSILON);
//        Assert.assertArrayEquals(expectedSp3Ds, Arrays.copyOfRange(esDsHi[2][1], 1885, 1895), EPSILON);
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
        esDsLoAtNode = q48One.getEsDsAtNode(nodeIdx, true);
        // printing
        // System.out.println("D's before prop in t: " + Arrays.toString(esDsHiAtNode[1])); // D's

        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i=0; i < esDsLoAtNode[ithDim].length; i++) {
                scratchAtNode[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }

        // just propagate in t, in place
        q48One.propagateTInPlace(esDsLoAtNode, scratchAtNode, true);

        esDsLoAtNode = q48One.getEsDsAtNode(nodeIdx, true);
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
        esDsHiAtNode = q48Two.getEsDsAtNode(nodeIdx, true);

        // making deep copy for assert below (checking that initial D's are correct)
        double[][] esDsHiAtNodeInitial = new double[esDsHiAtNode.length][esDsHiAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsHiAtNode[0], esDsHiAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsHiAtNode[1], esDsHiAtNode[1].length); // D

        // just propagate in x, in place
        double[] fftedfY = q48Two.getfY(true);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q48Two.getnXbins(true),0, 0, 1.0); // getting real part for assert below

        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsHiAtNode[0].length];
        q48Two.propagateXInPlace(esDsHiAtNode, scratchAtNode, true);

        esDsHiAtNode = q48Two.getEsDsAtNode(nodeIdx, true);
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

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 1024 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes1024BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsHiAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsHiAtNode = q1024.getEsDsAtNode(nodeIdx, true);

        String initialValuesStr = "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.07356011714376e-320, 1.07646593078779e-316, 2.19444210413816e-313, 4.29809786779161e-310, 8.0882896186969e-307, 1.46239706910102e-303, 2.54040022766109e-300, 4.24001310304921e-297, 6.79924162752444e-294, 1.04756739392715e-290, 1.55071393737006e-287, 2.20551303475434e-284, 3.01380943524079e-281, 3.95685529630673e-278, 4.99128960231452e-275, 6.04927865726296e-272, 7.0440532888614e-269, 7.88079255427205e-266, 8.47120866306929e-263, 8.74881199193118e-260, 8.68122493117834e-257, 8.27639395486764e-254, 7.58105280018574e-251, 6.67184746525239e-248, 5.64145225220522e-245, 4.58314867209543e-242, 3.5773809909959e-239, 2.68283933469824e-236, 1.93309125470737e-233, 1.33825312975345e-230, 8.90127872151017e-228, 5.68846058458759e-225, 3.49273251351755e-222, 2.06045968516333e-219, 1.16786094930117e-216, 6.35984267922332e-214, 3.32759014587752e-211, 1.67279032117111e-208, 8.0794426979051e-206, 3.74929423451527e-203, 1.67165196678454e-200, 7.16094641547923e-198, 2.9472922697571e-195, 1.16547834021394e-192, 4.42805931067358e-190, 1.61640755832608e-187, 5.66913157281163e-185, 1.91033890838957e-182, 6.18489778402691e-180, 1.92390158616881e-177, 5.74991466774607e-175, 1.65108019165656e-172, 4.55515495747242e-170, 1.20744223918252e-167, 3.07508973936726e-165, 7.52449650211712e-163, 1.76898974729278e-160, 3.99577851833656e-158, 8.67172946143426e-156, 1.8081658128976e-153, 3.62242170861341e-151, 6.97249132581543e-149, 1.28945199427945e-146, 2.2911345498459e-144, 3.91132795004428e-142, 6.41543470441782e-140, 1.01101166750794e-137, 1.53078594728388e-135, 2.22690139125134e-133, 3.11254563200686e-131, 4.17983067181575e-129, 5.39299317374918e-127, 6.68542888358834e-125, 7.96263661950535e-123, 9.11197964822431e-121, 1.00183639776346e-118, 1.0583007202168e-116, 1.07411207300394e-114, 1.04741388109478e-112, 9.8133044576004e-111, 8.83365515851868e-109, 7.64000830876149e-107, 6.34856310565053e-105, 5.06856793303167e-103, 3.88797374111874e-101, 2.86542862629871e-99, 2.02900953617604e-97, 1.38040588402555e-95, 9.02314083909614e-94, 5.66678703028705e-92, 3.41935559160492e-90, 1.98234784757337e-88, 1.10418967243195e-86, 5.90929564931781e-85, 3.03847716959186e-83, 1.50108213729065e-81, 7.12493910800289e-80, 3.249272073547e-78, 1.42370240784759e-76, 5.99350099634455e-75, 2.42420958981576e-73, 9.42080400617952e-72, 3.51749908518976e-70, 1.26185147112051e-68, 4.34921326859857e-67, 1.44026163054372e-65, 4.58247704739833e-64, 1.40083642686346e-62, 4.11436460605744e-61, 1.16103776130575e-59, 3.14787975955296e-58, 8.20008107166511e-57, 2.05232614558392e-55, 4.93517810313117e-54, 1.14021697818816e-52, 2.53104809320918e-51, 5.39810728877691e-50, 1.10614190996888e-48, 2.17775191065533e-47, 4.11940204481734e-46, 7.48666115977061e-45, 1.30728535506374e-43, 2.19321311877783e-42, 3.53524482050667e-41, 5.47502838470966e-40, 8.14669535505585e-39, 1.16467511994726e-37, 1.59976555140121e-36, 2.11123270049027e-35, 2.67697359850864e-34, 3.26122146967924e-33, 3.81719826927329e-32, 4.29276747132557e-31, 4.63829355451241e-30, 4.81512226367853e-29, 4.80269080001678e-28, 4.60246141769576e-27, 4.23763850701887e-26, 3.74874480468359e-25, 3.18622226540176e-24, 2.60192323984759e-23, 2.0414611188613e-22, 1.53891972534128e-21, 1.1146000045441e-20, 7.75622386349327e-20, 5.18572940220102e-19, 3.33117606475986e-18, 2.05595471433372e-17, 1.21915162591241e-16, 6.94592549713167e-16, 3.80216307581598e-15, 1.99967574969938e-14, 1.01045421670732e-13, 4.90571057139242e-13, 2.28831298036031e-12, 1.02555072735932e-11, 4.41597992627408e-11, 1.82694408167278e-10, 7.26192300358372e-10, 2.7733599883306e-09, 1.01762805632897e-08, 3.58756781592795e-08, 1.21517656996468e-07, 3.95463928124892e-07, 1.23652410003313e-06, 3.7147236891104e-06, 1.07220706893955e-05, 2.9734390294686e-05, 7.92259818206398e-05, 0.000202817041309727, 0.000498849425801082, 0.0011788613551308, 0.00267660451529767, 0.00583893851582903, 0.0122380386022749, 0.0246443833694605, 0.0476817640292964, 0.0886369682387583, 0.158309031659594, 0.271659384673714, 0.447890605896856, 0.709491856924619, 1.07981933026374, 1.57900316601789, 2.21841669358911, 2.99454931271487, 3.88372109966421, 4.83941449038288, 5.79383105522965, 6.66449205783597, 7.36540280606644, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606649, 6.66449205783597, 5.79383105522965, 4.83941449038288, 3.8837210996643, 2.99454931271495, 2.21841669358911, 1.57900316601789, 1.07981933026378, 0.709491856924646, 0.447890605896856, 0.271659384673714, 0.158309031659602, 0.088636968238763, 0.0476817640292964, 0.0246443833694605, 0.0122380386022757, 0.00583893851582942, 0.00267660451529767, 0.0011788613551308, 0.000498849425801082, 0.000202817041309743, 7.92259818206398e-05, 2.9734390294686e-05, 1.07220706893955e-05, 3.71472368911076e-06, 1.23652410003313e-06, 3.95463928124892e-07, 1.21517656996468e-07, 3.58756781592834e-08, 1.01762805632909e-08, 2.7733599883306e-09, 7.26192300358372e-10, 1.82694408167301e-10, 4.41597992627465e-11, 1.02555072735932e-11, 2.28831298036031e-12, 4.9057105713931e-13, 1.01045421670746e-13, 1.99967574969938e-14, 3.80216307581598e-15, 6.94592549713273e-16, 1.2191516259126e-16, 2.05595471433372e-17, 3.33117606475986e-18, 5.18572940220102e-19, 7.7562238634946e-20, 1.1146000045441e-20, 1.53891972534128e-21, 2.0414611188613e-22, 2.60192323984807e-23, 3.18622226540176e-24, 3.74874480468359e-25, 4.23763850701887e-26, 4.60246141769668e-27, 4.80269080001775e-28, 4.81512226367853e-29, 4.63829355451241e-30, 4.29276747132649e-31, 3.81719826927412e-32, 3.26122146967924e-33, 2.67697359850864e-34, 2.11123270049075e-35, 1.59976555140158e-36, 1.16467511994726e-37, 8.14669535505585e-39, 5.47502838471098e-40, 3.53524482050753e-41, 2.19321311877783e-42, 1.30728535506374e-43, 7.48666115977061e-45, 4.11940204481841e-46, 2.17775191065533e-47, 1.10614190996888e-48, 5.39810728877691e-50, 2.53104809320987e-51, 1.14021697818816e-52, 4.93517810313117e-54, 2.05232614558392e-55, 8.20008107166747e-57, 3.14787975955388e-58, 1.16103776130575e-59, 4.11436460605744e-61, 1.40083642686389e-62, 4.58247704739973e-64, 1.44026163054372e-65, 4.34921326859857e-67, 1.26185147112091e-68, 3.51749908519088e-70, 9.42080400617952e-72, 2.42420958981576e-73, 5.99350099634653e-75, 1.42370240784806e-76, 3.249272073547e-78, 7.12493910800289e-80, 1.50108213729117e-81, 3.03847716959291e-83, 5.90929564931781e-85, 1.10418967243195e-86, 1.98234784757337e-88, 3.41935559160616e-90, 5.66678703028705e-92, 9.02314083909614e-94, 1.38040588402555e-95, 2.0290095361768e-97, 2.86542862629871e-99, 3.88797374111874e-101, 5.06856793303245e-103, 6.34856310565152e-105, 7.64000830876451e-107, 8.83365515851868e-109, 9.81330445760197e-111, 1.04741388109495e-112, 1.07411207300438e-114, 1.0583007202168e-116, 1.00183639776362e-118, 9.11197964822584e-121, 7.96263661950872e-123, 6.68542888358834e-125, 5.39299317375011e-127, 4.17983067181647e-129, 3.11254563200822e-131, 2.22690139125134e-133, 1.53078594728388e-135, 1.01101166750812e-137, 6.41543470441898e-140, 3.91132795004428e-142, 2.2911345498459e-144, 1.28945199427969e-146, 6.97249132581673e-149, 3.62242170861341e-151, 1.8081658128976e-153, 8.67172946143591e-156, 3.99577851833733e-158, 1.76898974729363e-160, 7.52449650211712e-163, 3.07508973936786e-165, 1.20744223918276e-167, 4.55515495747468e-170, 1.65108019165656e-172, 5.74991466774723e-175, 1.9239015861692e-177, 6.18489778403008e-180, 1.91033890838957e-182, 5.66913157281281e-185, 1.61640755832642e-187, 4.4280593106759e-190, 1.16547834021394e-192, 2.9472922697571e-195, 7.16094641548077e-198, 1.67165196678491e-200, 3.74929423451527e-203, 8.0794426979051e-206, 1.67279032117148e-208, 3.32759014587826e-211, 6.35984267922332e-214, 1.16786094930117e-216, 2.0604596851638e-219, 3.49273251351835e-222, 5.68846058459019e-225, 8.90127872151017e-228, 1.33825312975376e-230, 1.93309125470782e-233, 2.68283933470013e-236, 3.5773809909959e-239, 4.58314867209652e-242, 5.64145225220657e-245, 6.67184746525559e-248, 7.58105280018574e-251, 8.27639395486965e-254, 8.68122493118046e-257, 8.74881199193764e-260, 8.47120866306929e-263, 7.88079255427205e-266, 7.04405328886316e-269, 6.04927865726448e-272, 4.99128960231452e-275, 3.95685529630673e-278, 3.01380943524156e-281, 2.20551303475491e-284, 1.55071393737006e-287, 1.04756739392715e-290, 6.79924162752622e-294, 4.24001310305033e-297, 2.54040022766244e-300, 1.46239706910102e-303, 8.08828961869907e-307, 4.2980978677928e-310, 2.19444210413816e-313, 1.07646593078779e-316, 5.07356011714376e-320, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0";
        String[] initialValuesStrArray = initialValuesStr.split(", ");
        double[] initialValues = new double[initialValuesStrArray.length];
        double[] expectedInitialDs = new double[initialValuesStrArray.length];
        for (int i=0; i<initialValues.length; i++) {
            double rValue = Double.parseDouble(initialValuesStrArray[i]);
            // if (i>313) System.out.println("i = " + i + " Java = " + esDsHiAtNode[1][i] + " R = " + rValue);
            q1024.setEsDsAtNodeElementAtDim(0, 1, i, rValue, true); // D's to fix rounding differences and match R and PropagatesQuaSSETest.testPropagateChOneCh1024QuaSSETest
            expectedInitialDs[i] = Double.parseDouble(initialValuesStrArray[i]); // D's for assert
        } // getting initial values

        // making deep copy for assert below (checking that initial D's are correct)
        double[][] esDsHiAtNodeInitial = new double[esDsHiAtNode.length][esDsHiAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsHiAtNode[0], esDsHiAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsHiAtNode[1], esDsHiAtNode[1].length); // D

        // just propagate in x, in place
        double[] fftedfY = q1024.getfY(true);
        double[] realFFTedfY = new double[fftedfY.length]; // just for test, not used in propagate in X
        everyOtherInPlace(fftedfY, realFFTedfY, q1024.getnXbins(true),0, 0, 1.0); // getting real part for assert below

        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = q1024.getScratchAtNode(0, true);
        q1024.propagateXInPlace(esDsHiAtNode, scratchAtNode, true);

        esDsHiAtNode = q1024.getEsDsAtNode(nodeIdx, true);
        double[] esLoAtNode = esDsHiAtNode[0];
        double[] dsLoAtNode = esDsHiAtNode[1];

        System.out.println("Final D's: " + Arrays.toString(dsLoAtNode));

        double[] expectedFFTedfY = new double[] { 1, 0.999999749692906, 0.999998998781049, 0.9999977472927, 0.999995995274977, 0.999993742793842, 0.999990989934102, 0.999987736799399, 0.999983983512212, 0.999979730213853, 0.999974977064454, 0.999969724242971, 0.99996397194717, 0.999957720393623, 0.999950969817696, 0.999943720473548, 0.999935972634113, 0.999927726591093, 0.999918982654949, 0.999909741154885, 0.999900002438841, 0.999889766873476, 0.999879034844152, 0.999867806754928, 0.999856083028535, 0.999843864106368, 0.999831150448462, 0.999817942533483, 0.999804240858701, 0.99979004593998, 0.999775358311752, 0.999760178527001, 0.999744507157237, 0.999728344792482, 0.999711692041242, 0.999694549530485, 0.999676917905621, 0.999658797830471, 0.999640189987249, 0.999621095076532, 0.999601513817236, 0.999581446946586, 0.999560895220092, 0.999539859411516, 0.999518340312849, 0.999496338734275, 0.999473855504143, 0.999450891468938 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsFirst48AfterPropX = new double[] { 0, 0, 0, 0, 1.11022302462516e-16, 0, 0, 0, 2.77555756156289e-17, 2.77555756156289e-17, 0, 8.32667268468867e-17, 0, 1.04083408558608e-16, 5.89805981832114e-17, 0, 0, 4.99275100429575e-17, 1.74936020530536e-16, 0, 5.47996435555642e-17, 1.64060117487791e-16, 2.55997708857483e-16, 0, 8.5609964086127e-17, 9.49329645689794e-17, 2.91487168428854e-16, 3.34092533592135e-17, 0, 2.19845474956207e-16, 3.57227590802606e-16, 3.37917798902023e-16, 0, 4.14102761209603e-17, 0, 0, 0, 0, 0, 0, 0, 3.49655200626575e-17, 2.86229373536173e-17, 0, 0, 0, 1.55257751099924e-16, 0 };
        // double[] expectedSp1DsLast48AfterPropX = new double[] {  };

        Assert.assertArrayEquals(Arrays.copyOfRange(expectedInitialDs, 315, 325), Arrays.copyOfRange(esDsHiAtNodeInitial[1], 315, 325), 1E-16);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(realFFTedfY, 0, 48), 1E-15);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), 1E-16);
        // Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), EPSILON2);
    }

    // now we'll test the integration inside the class
    // double[][][] esDsHi = q2.getEsDs(false);

    // TODO: in another test, do inside class
    // q2.processBranch(myTree2.getNode(nodeIdx));
}
