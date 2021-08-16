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

    static QuaSSEDistribution q, q2, q3, q4;
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
        data = Arrays.asList(-0.19537143, 0.00433218, 0.25996570);
        quTraitrp = new RealParameter();
        quTraitrp.initByName("value", data, "keys", spNames);

        String spNames2 = "sp1 sp2";
        data2 = Arrays.asList(0.0, 0.0);
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
//        q = new QuaSSEDistribution();
//        q.initByName("dt", dtrp, "tc", tcrp,
//                "nX", nXbinsip, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
//                "drift", driftrp, "diffusion", diffusionrp,
//                "q2mLambda", lfn, "q2mMu", cfn,
//                "tree", myTree,
//                "q2d", nfn);
//
        q2 = new QuaSSEDistribution();
        q2.initByName("dt", dtrp, "tc", tcrp,
                "nX", nXbinsip2, "dX", dxBinrp, "xMid", xMidrp, "flankWidthScaler", flankWidthScalerrp, "hiLoRatio", hiLoRatiorp,
                "drift", driftrp, "diffusion", diffusionrp,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", myTree2,
                "q2d", nfn2);

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
//    @Test
//    public void testInitializationOfTips() {
//        double[][][] esDsHi = q.getEsDs(false);
//
//        double[] expectedSp1Ds = new double[] { 1.60021460660503, 1.74807698146288, 1.90483403228819, 2.07046547520583, 2.24487988376146, 2.42790945985862, 2.61930535022212, 2.81873362725308, 3.02577205091375, 3.23990772345705 };
//        double[] expectedSp2Ds = new double[] { 4.27279148709871e-07, 5.69933520491392e-07, 7.58317302503376e-07, 1.00644951269517e-06, 1.33243884785178e-06, 1.75961170550665e-06, 2.3179318461188e-06, 3.04578151100852e-06, 3.99218899881498e-06, 5.21960579763602e-06 };
//        double[] expectedSp3Ds = new double[] { 1.27621312438669e-25, 2.19814297381616e-25, 3.77661688851761e-25, 6.47238270601096e-25, 1.1064701711749e-24, 1.88681576360346e-24, 3.20947166037177e-24, 5.44567674756272e-24, 9.21689065275829e-24, 1.55607768879411e-23 };
//
//        Assert.assertArrayEquals(expectedSp1Ds, Arrays.copyOfRange(esDsHi[0][1], 1885, 1895), EPSILON);
//        Assert.assertArrayEquals(expectedSp2Ds, Arrays.copyOfRange(esDsHi[1][1], 1885, 1895), EPSILON);
//        Assert.assertArrayEquals(expectedSp3Ds, Arrays.copyOfRange(esDsHi[2][1], 1885, 1895), EPSILON);
//    }

    /*
     * Checks that propagate methods in time for E's and D's
     * inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 1/20 = 0.05
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0
     */
    @Test
    public void testIntegrateOneBranchHiResOutsideClassJustT() {

        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsHiAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsHiAtNode = q2.getEsDsAtNode(nodeIdx, true);
        // printing
        System.out.println("D's before prop in t: " + Arrays.toString(esDsHiAtNode[1])); // D's

        double[][] scratchAtNode = new double[2][esDsHiAtNode[0].length];
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i=0; i < esDsHiAtNode[ithDim].length; i++) {
                scratchAtNode[ithDim][i] = esDsHiAtNode[ithDim][i];
            }
        }

        // just propagate in t, in place
        q2.propagateTInPlace(esDsHiAtNode, scratchAtNode, true);

        esDsHiAtNode = q2.getEsDsAtNode(nodeIdx, true);
        double[] esHiAtNode = esDsHiAtNode[0];
        double[] dsHiAtNode = esDsHiAtNode[1];

        // printing
         for (int i=0; i<esDsHiAtNode[0].length; ++i) {
             System.out.println("e" + i + " = " + esHiAtNode[i] + " d" + i + " = " + dsHiAtNode[i]);
         }

        double[] expectedSp1EsAfterPropT = new double[] { 0.00149370994379142, 0.00149368786296715, 0.00149366566009926, 0.00149364334116711, 0.0014936209122861, 0.00149359837970088, 0.00149357574977893, 0.00149355302900443, 0.0014935302239704, 0.00149350734137188, 0.00149348438799768, 0.00149346137072356, 0.00149343829650327, 0.00149341517236111, 0.00149339200538258, 0.0014933688027066, 0.00149334557151628, 0.0014933223190299, 0.00149329905249347, 0.00149327577916897, 0.00149325250632815, 0.00149322924124165, 0.00149320599116985, 0.00149318276335575, 0.00149315956501302, 0.00149313640331981, 0.00149311328540833, 0.00149309021835688, 0.00149306720918005, 0.00149304426482232, 0.00149302139214795, 0.00149299859793396, 0.00149297588886225, 0.00149295327151184, 0.00149293075235137, 0.00149290833773301, 0.00149288603388406, 0.00149286384690245, 0.00149284178274917, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropT = new double[] { 0.00579005736174394, 0.0121352277429472, 0.0244366205695344, 0.0472783719250338, 0.0878844494649893, 0.156960263942535, 0.269336720628206, 0.444047657299951, 0.703382841402857, 1.07048880941099, 1.5653111283757, 2.19911224723919, 2.96839922078803, 3.84968692464613, 4.79685632357302, 5.7427010526194, 6.60547283589975, 7.29994891700679, 7.75111077151784, 7.90744686571959, 7.75062709006449, 7.29903803124074, 6.60423685237185, 5.74126889110111, 4.79536172096607, 3.84824841117278, 2.96710606450487, 2.19801824876123, 1.56443588948085, 1.06982441610763, 0.702903170869044, 0.443717712813496, 0.269120201314488, 0.156824570265688, 0.0878031690333764, 0.0472318065561379, 0.0244110920425363, 0.0121218288516093, 0.00578332208315072, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esHiAtNode, 0, 48), EPSILON2);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(dsHiAtNode, 0, 48), EPSILON2);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value  for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 1/20 = 0.05
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0
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
