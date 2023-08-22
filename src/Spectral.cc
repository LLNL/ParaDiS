/*---------------------------------------------------------------------------
 *
 *  This module implements functions to perform stress / force calculations
 *  using the DDD spectral method (DDD-FFT).
 *  See Bertin et al., MSMSE (2015); Bertin, IJP (2019)
 *
 *  Author: Nicolas Bertin
 *  bertin1@llnl.gov
 *
 *-------------------------------------------------------------------------*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Decomp.h"
#include "Spectral.h"

#define NS_DISTRIB         1
#define SINGLY_CONVOLUTED  0
#define DISP_GRAD          0
#define MAXCELL2PERCELL   20

#ifdef SPECTRAL

/*-------------------------------------------------------------------------
 *
 *      Function:     GaussQuadCoeffAll
 *
 *------------------------------------------------------------------------*/
void GaussQuadCoeffAll(int intOrder, real8 *positions, real8 *weights)
{
        int i;

        if (intOrder == 1) {
            positions[0] = 0.0;
            weights[0] = 2.0;
        } else if (intOrder == 2) {
            positions[0] = -0.577350269189626;
            positions[1] = -positions[0];
            weights[0] = 1.0;
            weights[1] = 1.0;
        } else if (intOrder == 3) {
            positions[0] = -0.774596669241483;
            positions[1] = 0.0;
            positions[2] = -positions[0];
            weights[0] = 5.0/9.0;
            weights[1] = 8.0/9.0;
            weights[2] = weights[0];
        } else if (intOrder == 4) {
            positions[0] = -0.861136311594053;
            positions[1] = -0.339981043584856;
            positions[2] = -positions[1];
            positions[3] = -positions[0];
            weights[0] = 0.347854845137454;
            weights[1] = 0.652145154862546;
            weights[2] = weights[1];
            weights[3] = weights[0];
        } else if (intOrder == 5) {
            positions[0] = -0.906179845938664;
            positions[1] = -0.538469310105683;
            positions[2] = 0.0;
            positions[3] = -positions[1];
            positions[4] = -positions[0];
            weights[0] = 0.236926885056189;
            weights[1] = 0.478628670499366;
            weights[2] = 0.568888888888888;
            weights[3] = weights[1];
            weights[4] = weights[0];
        } else if (intOrder == 6) {
            positions[0] = -0.932469514203152;
            positions[1] = -0.661209386466265;
            positions[2] = -0.238619186083197;
            positions[3] = -positions[2];
            positions[4] = -positions[1];
            positions[5] = -positions[0];
            weights[0] = 0.171324492379170;
            weights[1] = 0.360761573048139;
            weights[2] = 0.467913934572691;
            weights[3] = weights[2];
            weights[4] = weights[1];
            weights[5] = weights[0];
		} else if (intOrder == 8) {
			positions[0] = -0.960289856497536;	weights[0] = 	0.101228536290376;
			positions[1] = -0.796666477413627;	weights[1] = 	0.222381034453374;
			positions[2] = -0.525532409916329;	weights[2] = 	0.313706645877887;
			positions[3] = -0.183434642495650;	weights[3] = 	0.362683783378362;
			positions[4] =  0.183434642495650;	weights[4] = 	0.362683783378362;
			positions[5] =  0.525532409916329;	weights[5] = 	0.313706645877887;
			positions[6] =  0.796666477413627;	weights[6] = 	0.222381034453374;
			positions[7] =  0.960289856497536;	weights[7] = 	0.101228536290376;
		} else if (intOrder == 16) {
			positions[0] = -0.989400934991650;	weights[0] = 	0.027152459411754;
			positions[1] = -0.944575023073233;	weights[1] = 	0.062253523938648;
			positions[2] = -0.865631202387832;	weights[2] = 	0.095158511682493;
			positions[3] = -0.755404408355003;	weights[3] = 	0.124628971255534;
			positions[4] = -0.617876244402644;	weights[4] = 	0.149595988816577;
			positions[5] = -0.458016777657227;	weights[5] = 	0.169156519395003;
			positions[6] = -0.281603550779259;	weights[6] = 	0.182603415044924;
			positions[7] = -0.095012509837637;	weights[7] = 	0.189450610455068;
			positions[8] =  0.095012509837637;	weights[8] = 	0.189450610455068;
			positions[9] = 0.281603550779259;	weights[9] = 	0.182603415044924;
			positions[10] = 0.458016777657227;	weights[10] = 	0.169156519395003;
			positions[11] = 0.617876244402644;	weights[11] = 	0.149595988816577;
			positions[12] = 0.755404408355003;	weights[12] = 	0.124628971255534;
			positions[13] = 0.865631202387832;	weights[13] = 	0.095158511682493;
			positions[14] = 0.944575023073233;	weights[14] = 	0.062253523938648;
			positions[15] = 0.989400934991650;	weights[15] = 	0.027152459411754;
		} else if (intOrder == 32) {
			positions[0] =   -.997263861849482;    weights[0] =    .007018610009470;
			positions[1] =   -.985611511545269;    weights[1] =    .016274394730906;
			positions[2] =   -.964762255587506;    weights[2] =    .025392065309262;
			positions[3] =   -.934906075937740;    weights[3] =    .034273862913020;
			positions[4] =   -.896321155766052;    weights[4] =    .042835898022227;
			positions[5] =   -.849367613732570;    weights[5] =    .050998059262376;
			positions[6] =   -.794483795967942;    weights[6] =    .058684093478535;
			positions[7] =   -.732182118740290;    weights[7] =    .065822222776362;
			positions[8] =   -.663044266930215;    weights[8] =    .072345794108850;
			positions[9] =   -.587715757240762;    weights[9] =    .078193895787070;
			positions[10] =   -.506899908932229;    weights[10] =    .083311924226945;
			positions[11] =   -.421351276130636;    weights[11] =    .087652093004404;
			positions[12] =   -.331868602282128;    weights[12] =    .091173878695764;
			positions[13] =   -.239287362252137;    weights[13] =    .093844399080805;
			positions[14] =   -.144471961582796;    weights[14] =    .095638720079274;
			positions[15] =   -.048307665687738;    weights[15] =    .096540088514728;
			positions[16] =    .048307665687738;    weights[16] =    .096540088514728;
			positions[17] =    .144471961582797;    weights[17] =    .095638720079275;
			positions[18] =    .239287362252137;    weights[18] =    .093844399080805;
			positions[19] =    .331868602282128;    weights[19] =    .091173878695764;
			positions[20] =    .421351276130635;    weights[20] =    .087652093004403;
			positions[21] =    .506899908932229;    weights[21] =    .083311924226947;
			positions[22] =    .587715757240762;    weights[22] =    .078193895787070;
			positions[23] =    .663044266930215;    weights[23] =    .072345794108848;
			positions[24] =    .732182118740289;    weights[24] =    .065822222776362;
			positions[25] =    .794483795967942;    weights[25] =    .058684093478536;
			positions[26] =    .849367613732570;    weights[26] =    .050998059262375;
			positions[27] =    .896321155766052;    weights[27] =    .042835898022227;
			positions[28] =    .934906075937740;    weights[28] =    .034273862913022;
			positions[29] =    .964762255587506;    weights[29] =    .025392065309262;
			positions[30] =    .985611511545268;    weights[30] =    .016274394730905;
			positions[31] =    .997263861849481;    weights[31] =    .007018610009470;
		} else if (intOrder == 64) {
			positions[0] =   -.999305041735772;    weights[0] =    .001783280721697;
			positions[1] =   -.996340116771955;    weights[1] =    .004147033260562;
			positions[2] =   -.991013371476744;    weights[2] =    .006504457968979;
			positions[3] =   -.983336253884626;    weights[3] =    .008846759826364;
			positions[4] =   -.973326827789911;    weights[4] =    .011168139460131;
			positions[5] =   -.961008799652054;    weights[5] =    .013463047896719;
			positions[6] =   -.946411374858403;    weights[6] =    .015726030476023;
			positions[7] =   -.929569172131939;    weights[7] =    .017951715775697;
			positions[8] =   -.910522137078502;    weights[8] =    .020134823153531;
			positions[9] =   -.889315445995114;    weights[9] =    .022270173808383;
			positions[10] =   -.865999398154093;    weights[10] =    .024352702568711;
			positions[11] =   -.840629296252580;    weights[11] =    .026377469715055;
			positions[12] =   -.813265315122797;    weights[12] =    .028339672614259;
			positions[13] =   -.783972358943341;    weights[13] =    .030234657072403;
			positions[14] =   -.752819907260532;    weights[14] =    .032057928354852;
			positions[15] =   -.719881850171611;    weights[15] =    .033805161837141;
			positions[16] =   -.685236313054233;    weights[16] =    .035472213256883;
			positions[17] =   -.648965471254657;    weights[17] =    .037055128540240;
			positions[18] =   -.611155355172394;    weights[18] =    .038550153178615;
			positions[19] =   -.571895646202634;    weights[19] =    .039953741132720;
			positions[20] =   -.531279464019895;    weights[20] =    .041262563242625;
			positions[21] =   -.489403145707052;    weights[21] =    .042473515123654;
			positions[22] =   -.446366017253464;    weights[22] =    .043583724529323;
			positions[23] =   -.402270157963991;    weights[23] =    .044590558163757;
			positions[24] =   -.357220158337668;    weights[24] =    .045491627927418;
			positions[25] =   -.311322871990211;    weights[25] =    .046284796581314;
			positions[26] =   -.264687162208767;    weights[26] =    .046968182816209;
			positions[27] =   -.217423643740007;    weights[27] =    .047540165714831;
			positions[28] =   -.169644420423993;    weights[28] =    .047999388596458;
			positions[29] =   -.121462819296121;    weights[29] =    .048344762234804;
			positions[30] =   -.072993121787799;    weights[30] =    .048575467441504;
			positions[31] =   -.024350292663424;    weights[31] =    .048690957009140;
			positions[32] =    .024350292663425;    weights[32] =    .048690957009139;
			positions[33] =    .072993121787799;    weights[33] =    .048575467441503;
			positions[34] =    .121462819296120;    weights[34] =    .048344762234803;
			positions[35] =    .169644420423993;    weights[35] =    .047999388596458;
			positions[36] =    .217423643740007;    weights[36] =    .047540165714831;
			positions[37] =    .264687162208767;    weights[37] =    .046968182816210;
			positions[38] =    .311322871990211;    weights[38] =    .046284796581314;
			positions[39] =    .357220158337668;    weights[39] =    .045491627927417;
			positions[40] =    .402270157963991;    weights[40] =    .044590558163758;
			positions[41] =    .446366017253464;    weights[41] =    .043583724529324;
			positions[42] =    .489403145707053;    weights[42] =    .042473515123653;
			positions[43] =    .531279464019895;    weights[43] =    .041262563242624;
			positions[44] =    .571895646202634;    weights[44] =    .039953741132721;
			positions[45] =    .611155355172393;    weights[45] =    .038550153178616;
			positions[46] =    .648965471254657;    weights[46] =    .037055128540241;
			positions[47] =    .685236313054233;    weights[47] =    .035472213256881;
			positions[48] =    .719881850171611;    weights[48] =    .033805161837142;
			positions[49] =    .752819907260532;    weights[49] =    .032057928354852;
			positions[50] =    .783972358943342;    weights[50] =    .030234657072402;
			positions[51] =    .813265315122797;    weights[51] =    .028339672614259;
			positions[52] =    .840629296252580;    weights[52] =    .026377469715055;
			positions[53] =    .865999398154092;    weights[53] =    .024352702568711;
			positions[54] =    .889315445995114;    weights[54] =    .022270173808382;
			positions[55] =    .910522137078503;    weights[55] =    .020134823153530;
			positions[56] =    .929569172131940;    weights[56] =    .017951715775697;
			positions[57] =    .946411374858403;    weights[57] =    .015726030476025;
			positions[58] =    .961008799652054;    weights[58] =    .013463047896719;
			positions[59] =    .973326827789911;    weights[59] =    .011168139460130;
			positions[60] =    .983336253884626;    weights[60] =    .008846759826364;
			positions[61] =    .991013371476744;    weights[61] =    .006504457968979;
			positions[62] =    .996340116771955;    weights[62] =    .004147033260562;
			positions[63] =    .999305041735772;    weights[63] =    .001783280721697;
		} else if (intOrder == 128) {
			positions[0] =   -.999824887947132;    weights[0] =    .000449380960291;
			positions[1] =   -.999077459977377;    weights[1] =    .001045812679341;
			positions[2] =   -.997733248625515;    weights[2] =    .001642503018669;
			positions[3] =   -.995792758534981;    weights[3] =    .002238288430963;
			positions[4] =   -.993257112900213;    weights[4] =    .002832751471457;
			positions[5] =   -.990127818491735;    weights[5] =    .003425526040911;
			positions[6] =   -.986406742724587;    weights[6] =    .004016254983738;
			positions[7] =   -.982096108435719;    weights[7] =    .004604584256703;
			positions[8] =   -.977198491463908;    weights[8] =    .005190161832676;
			positions[9] =   -.971716818747136;    weights[9] =    .005772637542867;
			positions[10] =   -.965654366431966;    weights[10] =    .006351663161706;
			positions[11] =   -.959014757853700;    weights[11] =    .006926892566899;
			positions[12] =   -.951801961341265;    weights[12] =    .007497981925635;
			positions[13] =   -.944020287830221;    weights[13] =    .008064589890485;
			positions[14] =   -.935674388277917;    weights[14] =    .008626377798618;
			positions[15] =   -.926769250878948;    weights[15] =    .009183009871660;
			positions[16] =   -.917310198080961;    weights[16] =    .009734153415007;
			positions[17] =   -.907302883401757;    weights[17] =    .010279479015831;
			positions[18] =   -.896753288049158;    weights[18] =    .010818660739503;
			positions[19] =   -.885667717345397;    weights[19] =    .011351376324079;
			positions[20] =   -.874052796958032;    weights[20] =    .011877307372741;
			positions[21] =   -.861915468939548;    weights[21] =    .012396139543952;
			positions[22] =   -.849262987577969;    weights[22] =    .012907562739267;
			positions[23] =   -.836102915060907;    weights[23] =    .013411271288616;
			positions[24] =   -.822443116955644;    weights[24] =    .013906964132952;
			positions[25] =   -.808291757507913;    weights[25] =    .014394345004168;
			positions[26] =   -.793657294762193;    weights[26] =    .014873122602147;
			positions[27] =   -.778548475506411;    weights[27] =    .015343010768865;
			positions[28] =   -.762974330044094;    weights[28] =    .015803728659399;
			positions[29] =   -.746944166797062;    weights[29] =    .016255000909785;
			positions[30] =   -.730467566741909;    weights[30] =    .016696557801589;
			positions[31] =   -.713554377683587;    weights[31] =    .017128135423111;
			positions[32] =   -.696214708369514;    weights[32] =    .017549475827118;
			positions[33] =   -.678458922447719;    weights[33] =    .017960327185008;
			positions[34] =   -.660297632272646;    weights[34] =    .018360443937331;
			positions[35] =   -.641741692562308;    weights[35] =    .018749586940545;
			positions[36] =   -.622802193910585;    weights[36] =    .019127523609951;
			positions[37] =   -.603490456158549;    weights[37] =    .019494028058706;
			positions[38] =   -.583818021628764;    weights[38] =    .019848881232830;
			positions[39] =   -.563796648226618;    weights[39] =    .020191871042132;
			positions[40] =   -.543438302412811;    weights[40] =    .020522792486959;
			positions[41] =   -.522755152051176;    weights[41] =    .020841447780752;
			positions[42] =   -.501759559136145;    weights[42] =    .021147646468221;
			positions[43] =   -.480464072404172;    weights[43] =    .021441205539209;
			positions[44] =   -.458881419833553;    weights[44] =    .021721949538052;
			positions[45] =   -.437024501037105;    weights[45] =    .021989710668461;
			positions[46] =   -.414906379552274;    weights[46] =    .022244328893799;
			positions[47] =   -.392540275033268;    weights[47] =    .022485652032746;
			positions[48] =   -.369939555349859;    weights[48] =    .022713535850235;
			positions[49] =   -.347117728597636;    weights[49] =    .022927844143688;
			positions[50] =   -.324088435024413;    weights[50] =    .023128448824387;
			positions[51] =   -.300865438877677;    weights[51] =    .023315229994063;
			positions[52] =   -.277462620177904;    weights[52] =    .023488076016536;
			positions[53] =   -.253893966422694;    weights[53] =    .023646883584448;
			positions[54] =   -.230173564226660;    weights[54] =    .023791557781003;
			positions[55] =   -.206315590902079;    weights[55] =    .023922012136704;
			positions[56] =   -.182334305985337;    weights[56] =    .024038168681024;
			positions[57] =   -.158244042714225;    weights[57] =    .024139957989019;
			positions[58] =   -.134059199461188;    weights[58] =    .024227319222816;
			positions[59] =   -.109794231127643;    weights[59] =    .024300200167972;
			positions[60] =   -.085463640504515;    weights[60] =    .024358557264690;
			positions[61] =   -.061081969604139;    weights[61] =    .024402355633849;
			positions[62] =   -.036663790968733;    weights[62] =    .024431569097850;
			positions[63] =   -.012223698960616;    weights[63] =    .024446180196261;
			positions[64] =    .012223698960616;    weights[64] =    .024446180196263;
			positions[65] =    .036663790968734;    weights[65] =    .024431569097849;
			positions[66] =    .061081969604140;    weights[66] =    .024402355633850;
			positions[67] =    .085463640504516;    weights[67] =    .024358557264692;
			positions[68] =    .109794231127644;    weights[68] =    .024300200167971;
			positions[69] =    .134059199461188;    weights[69] =    .024227319222816;
			positions[70] =    .158244042714225;    weights[70] =    .024139957989019;
			positions[71] =    .182334305985337;    weights[71] =    .024038168681024;
			positions[72] =    .206315590902079;    weights[72] =    .023922012136703;
			positions[73] =    .230173564226660;    weights[73] =    .023791557781003;
			positions[74] =    .253893966422694;    weights[74] =    .023646883584448;
			positions[75] =    .277462620177904;    weights[75] =    .023488076016536;
			positions[76] =    .300865438877677;    weights[76] =    .023315229994063;
			positions[77] =    .324088435024414;    weights[77] =    .023128448824387;
			positions[78] =    .347117728597636;    weights[78] =    .022927844143687;
			positions[79] =    .369939555349859;    weights[79] =    .022713535850236;
			positions[80] =    .392540275033268;    weights[80] =    .022485652032746;
			positions[81] =    .414906379552275;    weights[81] =    .022244328893799;
			positions[82] =    .437024501037104;    weights[82] =    .021989710668461;
			positions[83] =    .458881419833552;    weights[83] =    .021721949538052;
			positions[84] =    .480464072404172;    weights[84] =    .021441205539208;
			positions[85] =    .501759559136144;    weights[85] =    .021147646468223;
			positions[86] =    .522755152051176;    weights[86] =    .020841447780751;
			positions[87] =    .543438302412810;    weights[87] =    .020522792486962;
			positions[88] =    .563796648226618;    weights[88] =    .020191871042131;
			positions[89] =    .583818021628763;    weights[89] =    .019848881232830;
			positions[90] =    .603490456158548;    weights[90] =    .019494028058707;
			positions[91] =    .622802193910585;    weights[91] =    .019127523609950;
			positions[92] =    .641741692562307;    weights[92] =    .018749586940545;
			positions[93] =    .660297632272646;    weights[93] =    .018360443937330;
			positions[94] =    .678458922447719;    weights[94] =    .017960327185009;
			positions[95] =    .696214708369514;    weights[95] =    .017549475827118;
			positions[96] =    .713554377683587;    weights[96] =    .017128135423111;
			positions[97] =    .730467566741908;    weights[97] =    .016696557801589;
			positions[98] =    .746944166797062;    weights[98] =    .016255000909785;
			positions[99] =    .762974330044094;    weights[99] =    .015803728659399;
			positions[100] =    .778548475506412;    weights[100] =    .015343010768867;
			positions[101] =    .793657294762193;    weights[101] =    .014873122602146;
			positions[102] =    .808291757507914;    weights[102] =    .014394345004167;
			positions[103] =    .822443116955644;    weights[103] =    .013906964132951;
			positions[104] =    .836102915060907;    weights[104] =    .013411271288616;
			positions[105] =    .849262987577969;    weights[105] =    .012907562739268;
			positions[106] =    .861915468939549;    weights[106] =    .012396139543951;
			positions[107] =    .874052796958032;    weights[107] =    .011877307372740;
			positions[108] =    .885667717345397;    weights[108] =    .011351376324080;
			positions[109] =    .896753288049158;    weights[109] =    .010818660739504;
			positions[110] =    .907302883401757;    weights[110] =    .010279479015832;
			positions[111] =    .917310198080960;    weights[111] =    .009734153415007;
			positions[112] =    .926769250878948;    weights[112] =    .009183009871662;
			positions[113] =    .935674388277916;    weights[113] =    .008626377798617;
			positions[114] =    .944020287830220;    weights[114] =    .008064589890486;
			positions[115] =    .951801961341264;    weights[115] =    .007497981925634;
			positions[116] =    .959014757853700;    weights[116] =    .006926892566899;
			positions[117] =    .965654366431965;    weights[117] =    .006351663161707;
			positions[118] =    .971716818747137;    weights[118] =    .005772637542865;
			positions[119] =    .977198491463907;    weights[119] =    .005190161832677;
			positions[120] =    .982096108435719;    weights[120] =    .004604584256702;
			positions[121] =    .986406742724586;    weights[121] =    .004016254983739;
			positions[122] =    .990127818491734;    weights[122] =    .003425526040910;
			positions[123] =    .993257112900213;    weights[123] =    .002832751471459;
			positions[124] =    .995792758534981;    weights[124] =    .002238288430962;
			positions[125] =    .997733248625514;    weights[125] =    .001642503018669;
			positions[126] =    .999077459977376;    weights[126] =    .001045812679341;
			positions[127] =    .999824887947132;    weights[127] =    .000449380960292;
        } else {
			Fatal("Gaussian coefficients not available for intOrder = %d", intOrder);
		}

        for (i = 0; i < intOrder; i++) weights[i] *= 0.5;

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     BoxIntersection
 * 		Description:  Test if segment intersects an axis-aligned box
 *
 *-------------------------------------------------------------------------*/
int BoxIntersection(double p1[3], double p2[3], double box[3][2])
{
		int     i;
		double  c[3], e[3], m[3], d[3], ad[3];
		double  eps = 1.e-8;

		for (i = 0; i < 3; i++) {
			c[i] = 0.5*(box[i][0] + box[i][1]);
			e[i] = box[i][1] - c[i];
			m[i] = 0.5*(p1[i] + p2[i]);
			d[i] = p1[i] - m[i];
			m[i] -= c[i];
			ad[i] = fabs(d[i]);

			if (fabs(m[i]) > e[i] + ad[i]) return 0;
		}

		if (fabs(m[1]*d[2] - m[2]*d[1]) > e[1]*ad[2] + e[2]*ad[1]) return 0;
		if (fabs(m[2]*d[0] - m[0]*d[2]) > e[0]*ad[2] + e[2]*ad[0]) return 0;
		if (fabs(m[0]*d[1] - m[1]*d[0]) > e[0]*ad[1] + e[1]*ad[0]) return 0;

		return 1;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     OutsideBox
 * 		Description:  Test if point is outside of an axis-aligned box
 *
 *-------------------------------------------------------------------------*/
int OutsideBox(double p[3], double box[3][2])
{
		if (p[0] < box[0][0] || p[0] > box[0][1]) return 1;
		if (p[1] < box[1][0] || p[1] > box[1][1]) return 1;
		if (p[2] < box[2][0] || p[2] > box[2][1]) return 1;

		return 0;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SegPlaneIntersection
 * 		Description:  Computes the intersection between a plane defined by
 *                    its normal and a point, and a segment defined by
 *                    two points
 *
 *-------------------------------------------------------------------------*/
int SegPlaneIntersection(double n[3], double p[3], double p1[3],
                         double p2[3], double pint[3])
{
		double  t[3], D, w[3], N, st;
		double  eps = 1.e-5;

		VECTOR_COPY(pint, p1);

		t[0] = p2[0] - p1[0];
		t[1] = p2[1] - p1[1];
		t[2] = p2[2] - p1[2];

		D = n[0]*t[0] + n[1]*t[1] + n[2]*t[2];
		if (fabs(D) < eps) return 0;

		w[0] = p1[0] - p[0];
		w[1] = p1[1] - p[1];
		w[2] = p1[2] - p[2];

		N = -n[0]*w[0] - n[1]*w[1] - n[2]*w[2];
		st = N/D;
		if (st < 0.0 || st > 1.0) return 0;

		pint[0] += st*t[0];
		pint[1] += st*t[1];
		pint[2] += st*t[2];
		return 1;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SegBoxIntersectionPosition
 * 		Description:  Determine the intersection position between segment
 *                    pin-pout and an axis-aligned box
 *
 *-------------------------------------------------------------------------*/
int SegBoxIntersectionPosition(double pin[3], double pout[3],
                               double box[3][2], double pint[3])
{
		int     i, j;
		double  p[3], n[3];
		double  eps = 1.e-10;

		VECTOR_COPY(pint, pout);

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 2; j++) {

				if (j == 0) {
					if (pout[i] >= box[i][0]) continue;
				} else {
					if (pout[i] <= box[i][1]) continue;
				}

				// Point on the surface
				VECTOR_ZERO(p);
				p[i] = box[i][j];

				// Normal to the surface
				VECTOR_ZERO(n);
				n[i] = 2.0*(j - 0.5);

				// Compute intersection
				if (!SegPlaneIntersection(n, p, pin, pout, pint)) continue;

				// Correct rounding errors
				if (fabs(p[i]-pint[i]) < eps) pint[i] = p[i];

				if (!OutsideBox(pint, box)) return 1;
				else VECTOR_COPY(pint, pout);
			}
		}

		return 0;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     LineBoxIntersectionPosition
 * 		Description:  Determine the intersection position between a line
 *                    and an axis-aligned box
 *
 *-------------------------------------------------------------------------*/
int LineBoxIntersectionPosition(double p[3], double d[3], double box[3][2],
                                double pint1[3], double pint2[3])
{
		int     i;
		double  tmin, tmax, dinv, t1, t2, tmp;
		double  eps = 1.e-10;

		VECTOR_COPY(pint1, p);
		VECTOR_COPY(pint2, p);

		tmin = 0.0;
		tmax = 1.e20;

		for (i = 0; i < 3; i++) {
			if (fabs(d[i]) < eps) {
				if (p[i] < box[i][0] || p[i] > box[i][1]) return 0;
			} else {
				dinv = 1.0/d[i];
				t1 = (box[i][0] - p[i])*dinv;
				t2 = (box[i][1] - p[i])*dinv;
				if (t1 > t2) {
					tmp = t2;
					t2 = t1;
					t1 = tmp;
				}
				if (t1 > tmin) tmin = t1;
				if (t2 < tmax) tmax = t2;
				if (tmin > tmax) return 0;
			}
		}

		pint1[0] += tmin*d[0];
		pint1[1] += tmin*d[1];
		pint1[2] += tmin*d[2];

		pint2[0] += tmax*d[0];
		pint2[1] += tmax*d[1];
		pint2[2] += tmax*d[2];

		return 1;
}

/*------------------------------------------------------------------------
 *
 *      Function:    BubbleSort
 *
 *-----------------------------------------------------------------------*/
void BubbleSort(int id[], double val[], int n)
{
	int i, j;
	for (i = 0; i < n-1; i++) {
		for (j = 0; j < n-i-1; j++) {
			if (val[j] > val[j+1]) {
				double vtmp = val[j];
				val[j] = val[j+1];
				val[j+1] = vtmp;
				int itmp = id[j];
				id[j] = id[j+1];
				id[j+1] = itmp;
			}
		}
	}
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FFTAssignSlice
 *      Description:  Assign the memory of the FFT slice spanning the domain
 *
 *-------------------------------------------------------------------------*/
void FFTAssignSlice(Home_t *home)
{
#ifdef PARALLEL
        int i, j, nc, d, p;
        Node_t *node, *nbr;
        double p1[3], p2[3];
        Param_t *param = home->param;
        Spectral_t *spectral = home->spectral;

        int numDomains = home->numDomains;
        int myDomain = home->myDomain;

        double Hx = param->Lx / spectral->N;
        double Dx = param->Lx / numDomains;

        double pxmin = param->minSideX + myDomain * Dx;
        double pxmax = pxmin + Dx;

/*
 *		Loop over all the local nodes to determine all remote domains
 *      intersected by segments owned by the current task.
 */
		for (i = 0; i < home->newNodeKeyPtr; i++) {

			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			nc = node->numNbrs;

			p1[0] = node->x;
			p1[1] = node->y;
			p1[2] = node->z;

            pxmin = MIN(p1[0], pxmin);
            pxmax = MAX(p1[0], pxmax);

			for (j = 0; j < nc; j++) {

				nbr = GetNeighborNode(home, node, j);
				if (nbr == (Node_t *)NULL) continue;

				/* Avoid double counting */
				if (OrderNodes(node, nbr) != -1) continue;

				p2[0] = nbr->x;
				p2[1] = nbr->y;
				p2[2] = nbr->z;

				PBCPOSITION(param, p1[0], p1[1], p1[2], &p2[0], &p2[1], &p2[2]);

                pxmin = MIN(p2[0], pxmin);
                pxmax = MAX(p2[0], pxmax);
            }
        }

        int Dmin = floor((pxmin - Hx - param->minSideX) / Dx);
        int Dmax = floor((pxmax + Hx - param->minSideX) / Dx);

        int Dspan = Dmax-Dmin+1;

        Dspan = MAX(Dmax-myDomain, myDomain-Dmin);
        if (Dspan >= numDomains/2)
            Fatal("Overlapping FFT domain communications");
        Dspan = 2*Dspan+1;

        //printf("Proc %d: Dspan = %d\n", myDomain, Dspan);

        MPI_Allreduce(MPI_IN_PLACE, &Dspan, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

/*
 *      Reallocate the local memory for the alpha / stress tensors. The size of
 *      the memory block must correspond to the span of dislocation segments owned
 *      by the current domain, and is therefore likely to exceed the default size
 *      of the local FFT slice.
 */
        int N = spectral->N;
        ptrdiff_t local_Nx = spectral->local_Nx;

        if (Dspan != spectral->Dspan) {
            spectral->Dspan = Dspan;

            for (i = 0; i < 9; i++) {
    			spectral->fields[i] = (fftcomplex_t*)realloc(spectral->fields[i],
                                      Dspan*local_Nx*N*N*sizeof(fftcomplex_t));
    			if (spectral->fields[i] == NULL)
                    Fatal("Error reallocating memory for the spectral fields array");

                spectral->alpha[i] = spectral->fields[i];
    		}

            int stressptr[6] = {0, 4, 8, 5, 2, 1};
            for (i = 0; i < 6; i++)
                spectral->stress[i] = (double*)spectral->fields[stressptr[i]];

#if DISP_GRAD
            if (param->localRotation) {
                for (i = 0; i < 9; i++) {
        			spectral->dispgrad[i] = (double*)realloc(spectral->dispgrad[i],
                                             Dspan*local_Nx*N*N*sizeof(double));
        			if (spectral->dispgrad[i] == NULL)
                        Fatal("malloc error for displacement gradient field array");
        		}
            }
#endif
        }
#endif
        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     AlphaTensor
 *
 *-------------------------------------------------------------------------*/
void AlphaTensor(Home_t *home)
{
		int     i, j, k, l, nc;
		int     kx, ky, kz, ib, jb, kb, N;
		int     imin, jmin, kmin, imax, jmax, kmax;
		int     in1, in2, is[3], i1, i2, start;
		double  H[3], Hmax;
		double  p1[3], p2[3], x1[3], x2[3], x0[3];
		double  b[3], t[3], alpha[9], dt;
		double  R[3], Rt, d[3], x10[3], x20[3];
		double  s1, s2, s[3], ss1, ss2;
		double  bc[3], bs[3], bst[3], dc, box[3][2];
		double  xmin, ymin, zmin, xmax, ymax, zmax;
		double  ltot, wtot, W, A, B, C, V;
		Node_t  *node, *nbr;
		Param_t *param;
		Spectral_t *spectral;

		param = home->param;
		spectral = home->spectral;
		N = spectral->N;

        int myDomain = home->myDomain;
#ifdef PARALLEL
        int numDomains = home->numDomains;
        int Dspan = spectral->Dspan;
        int Dmin = myDomain-Dspan/2;
        ptrdiff_t local_Nx = spectral->local_Nx;
#else
        int Dspan = 1;
        int Dmin = myDomain;
        int local_Nx = spectral->local_Nx;
#endif

		H[0] = param->Lx/N;
		H[1] = param->Ly/N;
		H[2] = param->Lz/N;
		Hmax = MAX(H[0], MAX(H[1], H[2]));
		V = H[0]*H[1]*H[2];
		//V = 1.0/N/N;
		//V = param->simVol;

		ltot = 0.0;
		wtot = 0.0;

		int nseg = 0;

		for (kx = 0; kx < Dspan*local_Nx; kx++)
			for (ky = 0; ky < N; ky++)
				for (kz = 0; kz < N; kz++)
					for (k = 0; k < 9; k++)
						spectral->alpha[k][kz+N*ky+N*N*kx] = fftcomplex_t(0.0);

/*
 *		Loop over all the local nodes
 */
		for (i = 0; i < home->newNodeKeyPtr; i++) {

			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			nc = node->numNbrs;

			p1[0] = node->x;
			p1[1] = node->y;
			p1[2] = node->z;

/*
 *			Loop over every segment attached to the node
 */
			for (j = 0; j < nc; j++) {

				nbr = GetNeighborNode(home, node, j);
				if (nbr == (Node_t *)NULL) continue;

				/* Avoid double counting */
				if (OrderNodes(node, nbr) != -1) continue;
                //if (NodeOwnsSeg(home, node, nbr) == 0) continue;

				nseg++;

				b[0] = node->burgX[j];
				b[1] = node->burgY[j];
				b[2] = node->burgZ[j];

				p2[0] = nbr->x;
				p2[1] = nbr->y;
				p2[2] = nbr->z;

				PBCPOSITION(param, p1[0], p1[1], p1[2], &p2[0], &p2[1], &p2[2]);

				t[0] = p2[0] - p1[0];
				t[1] = p2[1] - p1[1];
				t[2] = p2[2] - p1[2];
				dt = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
				if (dt < 1.e-10) continue;
				ltot += dt;

				t[0] /= dt;
				t[1] /= dt;
				t[2] /= dt;

				for (k = 0; k < 3; k++)
					for (l = 0; l < 3; l++)
						alpha[k*3+l] = b[k]*t[l]/V;


				/* Segment Bounding box */
				xmin = MIN(p1[0], p2[0]) - 0.5*H[0] - param->minSideX;
				xmax = MAX(p1[0], p2[0]) - 0.5*H[0] - param->minSideX;
				ymin = MIN(p1[1], p2[1]) - 0.5*H[1] - param->minSideY;
				ymax = MAX(p1[1], p2[1]) - 0.5*H[1] - param->minSideY;
				zmin = MIN(p1[2], p2[2]) - 0.5*H[2] - param->minSideZ;
				zmax = MAX(p1[2], p2[2]) - 0.5*H[2] - param->minSideZ;

				imin = floor(xmin/H[0]);
				imax = floor(xmax/H[0]) + 1;
				jmin = floor(ymin/H[1]);
				jmax = floor(ymax/H[1]) + 1;
				kmin = floor(zmin/H[2]);
				kmax = floor(zmax/H[2]) + 1;

/*
 *				Loop over regional boxes
 */
				for (ib = imin; ib <= imax; ib++) {
					for (jb = jmin; jb <= jmax; jb++) {
						for (kb = kmin; kb <= kmax; kb++) {

							/* Box center */
							double bc[3];
							bc[0] =  (ib + 0.5)*H[0] + param->minSideX;
							bc[1] =  (jb + 0.5)*H[1] + param->minSideY;
							bc[2] =  (kb + 0.5)*H[2] + param->minSideZ;

							/* Check that segment intersects box */
							bs[0] = bc[0] - p1[0];
							bs[1] = bc[1] - p1[1];
							bs[2] = bc[2] - p1[2];
							CrossVector(bs, t, bst);
							dc = Normal(bst);
							if (dc > sqrt(3.0)*Hmax) continue;

							box[0][0] = bc[0] - H[0];
							box[0][1] = bc[0] + H[0];
							box[1][0] = bc[1] - H[1];
							box[1][1] = bc[1] + H[1];
							box[2][0] = bc[2] - H[2];
							box[2][1] = bc[2] + H[2];
							if (!BoxIntersection(p1, p2, box)) continue;

							// Find which portion of the segment lies within the box
							in1 = !OutsideBox(p1, box);
							in2 = !OutsideBox(p2, box);
							if (in1 && in2) {
								VECTOR_COPY(x1, p1);
								VECTOR_COPY(x2, p2);
							} else if (in1) {
								VECTOR_COPY(x1, p1);
								SegBoxIntersectionPosition(p1, p2, box, x2);
							} else if (in2) {
								SegBoxIntersectionPosition(p2, p1, box, x1);
								VECTOR_COPY(x2, p2);
							} else {
								if (!LineBoxIntersectionPosition(p1, t, box, x1, x2)) {
									//Fatal("Cannot find line box intersection point");
                                    continue;
								}
							}

							d[0] = x2[0] - x1[0];
							d[1] = x2[1] - x1[1];
							d[2] = x2[2] - x1[2];
							if (Normal(d) < 1.e-10) continue;

							// Parametrize segment
							R[0] = bc[0] - x1[0];
							R[1] = bc[1] - x1[1];
							R[2] = bc[2] - x1[2];
							Rt = R[0]*t[0] + R[1]*t[1] + R[2]*t[2];
							x0[0] = x1[0] + Rt*t[0];
							x0[1] = x1[1] + Rt*t[1];
							x0[2] = x1[2] + Rt*t[2];
							d[0] = R[0] - Rt*t[0];
							d[1] = R[1] - Rt*t[1];
							d[2] = R[2] - Rt*t[2];
							x10[0] = x1[0] - x0[0];
							x10[1] = x1[1] - x0[1];
							x10[2] = x1[2] - x0[2];
							x20[0] = x2[0] - x0[0];
							x20[1] = x2[1] - x0[1];
							x20[2] = x2[2] - x0[2];
							s1 = x10[0]*t[0] + x10[1]*t[1] + x10[2]*t[2];
							s2 = x20[0]*t[0] + x20[1]*t[1] + x20[2]*t[2];

							if (fabs(t[0]) < 1.e-20) s[0] = s2+1; else s[0] = d[0]/t[0];
							if (fabs(t[1]) < 1.e-20) s[1] = s2+1; else s[1] = d[1]/t[1];
							if (fabs(t[2]) < 1.e-20) s[2] = s2+1; else s[2] = d[2]/t[2];

							// First term
							W = s2 - s1;

							// Second term Ai
							for (k = 0; k < 3; k++) {
								A = 0.0;
								if (s[k] > s1 && s[k] < s2) {
									A += fabs(0.5*(s[k]-s1)*(2*d[k]-t[k]*(s[k]+s1)));
									A += fabs(0.5*(s2-s[k])*(2*d[k]-t[k]*(s2+s[k])));
								} else {
									A += fabs(0.5*(s2-s1)*(2*d[k]-t[k]*(s2+s1)));
								}
								W -= 1.0/H[k]*A;
							}

							// Third term Bij
							for (k = 0; k < 3; k++) {
								i1 = k;
								if (i1 == 2) i2 = 0; else i2 = i1+1;
								int ss_len = 4;
								int id[] = {0, 1, 2, 3};
								double ss[] = {s1, s2, s[i1], s[i2]};
								BubbleSort(id, ss, ss_len);
								start = 0;
								B = 0.0;
								for (l = 0; l < ss_len-1; l++) {
									if (id[l] == 0) start = 1;
									if (start) {
										ss1 = ss[l];
										ss2 = ss[l+1];
										B += fabs(d[i1]*d[i2]*(ss2-ss1)
										     -0.5*(d[i1]*t[i2]+d[i2]*t[i1])*(ss2*ss2-ss1*ss1)
										     +1.0/3.0*t[i1]*t[i2]*(ss2*ss2*ss2-ss1*ss1*ss1));
									}
									if (id[l+1] == 1) break;
								}
								W += 1.0/H[i1]/H[i2]*B;
							}

							// Forth term Cijk
							int ss_len = 5;
							int id[] = {0, 1, 2, 3, 4};
							double ss[] = {s1, s2, s[0], s[1], s[2]};
							BubbleSort(id, ss, ss_len);
							start = 0;
							C = 0.0;
							for (l = 0; l < ss_len-1; l++) {
								if (id[l] == 0) start = 1;
								if (start) {
									ss1 = ss[l];
									ss2 = ss[l+1];
									C += fabs(d[0]*d[1]*d[2]*(ss2-ss1)
									     -0.5*(t[0]*d[1]*d[2]+d[0]*t[1]*d[2]+d[0]*d[1]*t[2])*(ss2*ss2-ss1*ss1)
									     +1.0/3.0*(d[0]*t[1]*t[2]+t[0]*d[1]*t[2]+t[0]*t[1]*d[2])*(ss2*ss2*ss2-ss1*ss1*ss1)
									     -0.25*(t[0]*t[1]*t[2])*(ss2*ss2*ss2*ss2-ss1*ss1*ss1*ss1));
								}
								if (id[l+1] == 1) break;
							}
							W -= 1.0/H[0]/H[1]/H[2]*C;
							wtot += W;

                            /* Box index */
                        #ifdef PARALLEL
                            kx = ib - Dmin*local_Nx;
                            if (kx < 0 || kx >= Dspan*local_Nx)
                                Fatal("Local box out of Dpsan on proc %d", myDomain);
                        #else
                            kx = ib % N;
							if (kx < 0) kx += N;
                        #endif
                            ky = jb % N;
							if (ky < 0) ky += N;
							kz = kb % N;
							if (kz < 0) kz += N;

							for (k = 0; k < 9; k++)
								spectral->alpha[k][kz+N*ky+N*N*kx] += fftcomplex_t(W*alpha[k]);

						}
					}
				}

			}
		}

		//printf("ltot = %e, wtot = %e, wtot/ltot = %e\n", ltot, wtot, wtot/ltot);
		//printf("dens = %e\n", ltot/param->simVol/param->burgMag/param->burgMag);

		//free(densityGrid);

#ifdef PARALLEL

        MPI_Request request;
        MPI_Status  status;

        fftcomplex_t *recvBuff = (fftcomplex_t*)malloc(N*N*local_Nx*sizeof(fftcomplex_t));

        for (int d = 0; d < Dspan; d++) {
            if (d == Dspan/2) continue;

            int send = (myDomain + Dspan/2 - d) % numDomains;
            if (send < 0) send += numDomains;

            int recv = (myDomain - Dspan/2 + d) % numDomains;
            if (recv < 0) recv += numDomains;

            for (int k = 0; k < 9; k++) {

                /* Send to MPI task send and receive from MPI task recv */
                MPI_Irecv((double*)recvBuff, 2*N*N*local_Nx, MPI_DOUBLE, recv, k, MPI_COMM_WORLD, &request);
                MPI_Send((double*)&spectral->alpha[k][(Dspan-d-1)*N*N*local_Nx], 2*N*N*local_Nx, MPI_DOUBLE, send, k, MPI_COMM_WORLD);
                MPI_Wait(&request, &status);

                /* Unpack buffer to increment local alpha tensor with remote contribution */
                for (kx = 0; kx < local_Nx; kx++)
        			for (ky = 0; ky < N; ky++)
        				for (kz = 0; kz < N; kz++)
                            spectral->alpha[k][kz+N*ky+N*N*((myDomain-Dmin)*local_Nx+kx)] += recvBuff[kz+N*ky+N*N*kx];
            }
        }

        free(recvBuff);
#endif
#if 0
        FILE *fp;
        char filename[128];

        sprintf(filename, "alpha_%d.dat", myDomain);
        fp = fopen(filename, "w");

        for (kx = 0; kx < local_Nx; kx++)
			for (ky = 0; ky < N; ky++)
				for (kz = 0; kz < N; kz++) {
					for (k = 0; k < 9; k++)
						fprintf(fp, "%e ", spectral->alpha[k][kz+N*ky+N*N*((myDomain-Dmin)*local_Nx+kx)].re);
                    fprintf(fp, "\n");
                }

        fclose(fp);
#endif
		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FFT3DForward
 *
 *-------------------------------------------------------------------------*/
void FFT3DForward(fftw_complex *in, fftw_complex *out, int nx, int ny, int nz)
{
		fftw_plan plan;
#ifdef PARALLEL
        plan = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
#else
		plan = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
		fftw_execute(plan);

		fftw_destroy_plan(plan);
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FFT3DBackward
 *
 *-------------------------------------------------------------------------*/
void FFT3DBackward(fftw_complex *in, fftw_complex *out, int nx, int ny, int nz)
{
		fftw_plan plan;
#ifdef PARALLEL
        plan = fftw_mpi_plan_dft_3d(nx, ny, nz, in, out, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
		plan = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
		fftw_execute(plan);

		fftw_destroy_plan(plan);
}

/*---------------------------------------------------------------------------
 *
 *      Function:     fftcomplex_t operators
 *
 *-------------------------------------------------------------------------*/
fftcomplex_t operator+(const fftcomplex_t a, const fftcomplex_t b) {
    return fftcomplex_t(a.re + b.re, a.im + b.im);
}

fftcomplex_t operator-(const fftcomplex_t a, const fftcomplex_t b) {
    return fftcomplex_t(a.re - b.re, a.im - b.im);
}

fftcomplex_t operator*(const fftcomplex_t a, const fftcomplex_t b)
{
    return fftcomplex_t(a.re*b.re-a.im*b.im, a.im*b.re+a.re*b.im);
}

fftcomplex_t operator*(const double c, const fftcomplex_t a)
{
    return fftcomplex_t(c*a.re, c*a.im);
}

fftcomplex_t operator/(const double c, const fftcomplex_t a)
{
    double d = a.re*a.re + a.im*a.im;
    return fftcomplex_t(c/d*a.re, -c/d*a.im);
}

fftcomplex_t conj(const fftcomplex_t a)
{
    return fftcomplex_t(a.re, -a.im);
}

/*---------------------------------------------------------------------------
 *
 *      Function:     InverseComplex33
 *
 *-------------------------------------------------------------------------*/
void InverseComplex33(fftcomplex_t M[3][3])
{
		fftcomplex_t a, b, c, d, e, f, g, h, i, com;

		a = M[0][0]; b = M[0][1]; c = M[0][2];
		d = M[1][0]; e = M[1][1]; f = M[1][2];
		g = M[2][0]; h = M[2][1]; i = M[2][2];

		com = 1.0/(a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g);

		M[0][0] = com*(e*i - f*h);
		M[0][1] = com*(c*h - b*i);
		M[0][2] = com*(b*f - c*e);
		M[1][0] = com*(f*g - d*i);
		M[1][1] = com*(a*i - c*g);
		M[1][2] = com*(c*d - a*f);
		M[2][0] = com*(d*h - e*g);
		M[2][1] = com*(b*g - a*h);
		M[2][2] = com*(a*e - b*d);

		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     StressFromAlpha
 *
 *-------------------------------------------------------------------------*/
void StressFromAlpha(Home_t *home)
{
		double  xk[3], nxk, nxk2, nxk2inv, H[3];
		Param_t *param;
		Spectral_t *spectral;

		param = home->param;
		spectral = home->spectral;
		int N = spectral->N;

		H[0] = param->Lx/N;
		H[1] = param->Ly/N;
		H[2] = param->Lz/N;

		int kmax = N/2 + N % 2;

        int myDomain = home->myDomain;
#ifdef PARALLEL
        int numDomains = home->numDomains;
        int Dspan = spectral->Dspan;
        int Dmin = myDomain-Dspan/2;
        ptrdiff_t local_Nx = spectral->local_Nx;
        ptrdiff_t local_kx_start = spectral->local_kx_start;
#else
        int Dspan = 1;
        int Dmin = myDomain;
        int local_Nx = spectral->local_Nx;
        int local_kx_start = spectral->local_kx_start;
#endif

		fftcomplex_t *Tau11k, *Tau12k, *Tau13k, *Tau21k, *Tau22k, *Tau23k, *Tau31k, *Tau32k, *Tau33k;

        Tau11k = &spectral->alpha[0][(myDomain-Dmin)*N*N*local_Nx];
		Tau12k = &spectral->alpha[1][(myDomain-Dmin)*N*N*local_Nx];
		Tau13k = &spectral->alpha[2][(myDomain-Dmin)*N*N*local_Nx];
		Tau21k = &spectral->alpha[3][(myDomain-Dmin)*N*N*local_Nx];
		Tau22k = &spectral->alpha[4][(myDomain-Dmin)*N*N*local_Nx];
		Tau23k = &spectral->alpha[5][(myDomain-Dmin)*N*N*local_Nx];
		Tau31k = &spectral->alpha[6][(myDomain-Dmin)*N*N*local_Nx];
		Tau32k = &spectral->alpha[7][(myDomain-Dmin)*N*N*local_Nx];
		Tau33k = &spectral->alpha[8][(myDomain-Dmin)*N*N*local_Nx];


		FFT3DForward((fftw_complex*)Tau11k, (fftw_complex*)Tau11k, N, N, N);
		FFT3DForward((fftw_complex*)Tau12k, (fftw_complex*)Tau12k, N, N, N);
		FFT3DForward((fftw_complex*)Tau13k, (fftw_complex*)Tau13k, N, N, N);
		FFT3DForward((fftw_complex*)Tau21k, (fftw_complex*)Tau21k, N, N, N);
		FFT3DForward((fftw_complex*)Tau22k, (fftw_complex*)Tau22k, N, N, N);
		FFT3DForward((fftw_complex*)Tau23k, (fftw_complex*)Tau23k, N, N, N);
		FFT3DForward((fftw_complex*)Tau31k, (fftw_complex*)Tau31k, N, N, N);
		FFT3DForward((fftw_complex*)Tau32k, (fftw_complex*)Tau32k, N, N, N);
		FFT3DForward((fftw_complex*)Tau33k, (fftw_complex*)Tau33k, N, N, N);


		for (int kx = 0; kx < local_Nx; kx++) {
			for (int ky = 0; ky < N; ky++) {
				for (int kz = 0; kz < N; kz++) {

                    int kxx = kx + local_kx_start;
                    int kyy = ky;
                    int kzz = kz;
					if (kxx >= kmax) kxx -= N;
					if (kyy >= kmax) kyy -= N;
					if (kzz >= kmax) kzz -= N;
					int kind = kz+ky*N+kx*N*N;

					if (kxx == 0 && kyy == 0 && kzz == 0) {
						Tau11k[kind] = fftcomplex_t(0.0, 0.0);
						Tau12k[kind] = fftcomplex_t(0.0, 0.0);
						Tau13k[kind] = fftcomplex_t(0.0, 0.0);
						Tau22k[kind] = fftcomplex_t(0.0, 0.0);
						Tau23k[kind] = fftcomplex_t(0.0, 0.0);
						Tau33k[kind] = fftcomplex_t(0.0, 0.0);
						continue;
					}

#if NS_DISTRIB
#if SINGLY_CONVOLUTED
					fftcomplex_t swk = csqrt(spectral->wk[kind]); // pre-compute sqrt I guess...
#else
					fftcomplex_t swk = spectral->wk[kind]; // apply doubly-convoluted kernel
#endif
					Tau11k[kind] *= swk;
					Tau12k[kind] *= swk;
					Tau13k[kind] *= swk;
					Tau21k[kind] *= swk;
					Tau22k[kind] *= swk;
					Tau23k[kind] *= swk;
					Tau31k[kind] *= swk;
					Tau32k[kind] *= swk;
					Tau33k[kind] *= swk;
#endif

					xk[0] = 2.0*M_PI*kxx/H[0]/N;
					xk[1] = 2.0*M_PI*kyy/H[1]/N;
					xk[2] = 2.0*M_PI*kzz/H[2]/N;

					nxk2 = xk[0]*xk[0] + xk[1]*xk[1] + xk[2]*xk[2];
					nxk2inv = 1.0/nxk2;

					fftcomplex_t cxk[3];
					cxk[0] = fftcomplex_t(0.0, xk[0]);
					cxk[1] = fftcomplex_t(0.0, xk[1]);
					cxk[2] = fftcomplex_t(0.0, xk[2]);

					fftcomplex_t K[3][3];
					K[0][0] = -nxk2inv*(Tau12k[kind]*cxk[2]-Tau13k[kind]*cxk[1]);
					K[0][1] = -nxk2inv*(Tau13k[kind]*cxk[0]-Tau11k[kind]*cxk[2]);
					K[0][2] = -nxk2inv*(Tau11k[kind]*cxk[1]-Tau12k[kind]*cxk[0]);

					K[1][0] = -nxk2inv*(Tau22k[kind]*cxk[2]-Tau23k[kind]*cxk[1]);
					K[1][1] = -nxk2inv*(Tau23k[kind]*cxk[0]-Tau21k[kind]*cxk[2]);
					K[1][2] = -nxk2inv*(Tau21k[kind]*cxk[1]-Tau22k[kind]*cxk[0]);

					K[2][0] = -nxk2inv*(Tau32k[kind]*cxk[2]-Tau33k[kind]*cxk[1]);
					K[2][1] = -nxk2inv*(Tau33k[kind]*cxk[0]-Tau31k[kind]*cxk[2]);
					K[2][2] = -nxk2inv*(Tau31k[kind]*cxk[1]-Tau32k[kind]*cxk[0]);

					fftcomplex_t T[3][3];
					for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
						T[i][j] = fftcomplex_t(0.0, 0.0);
						for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
							T[i][j] += spectral->C[i][j][k][l]*K[k][l];
						}
					}

					fftcomplex_t ac[3][3];
					for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) {
						ac[i][k] = fftcomplex_t(0.0, 0.0);
						for (int j = 0; j < 3; j++) for (int l = 0; l < 3; l++) {
							ac[i][k] += spectral->C[k][j][i][l]*cxk[j]*conj(cxk[l]);
						}
					}
					InverseComplex33(ac);

					fftcomplex_t G[3][3][3][3];
					for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
						for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
							G[i][j][k][l] = ac[i][k]*cxk[j]*conj(cxk[l]);
						}
					}

					fftcomplex_t U[3][3];
					for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
						U[i][j] = fftcomplex_t(0.0, 0.0);
						for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
							U[i][j] -= G[i][j][k][l]*T[k][l];
						}
					}

#if DISP_GRAD
                    if (param->localRotation) {
    					Tau11k[kind] = U[0][0] + K[0][0];
    					Tau12k[kind] = U[0][1] + K[0][1];
    					Tau13k[kind] = U[0][2] + K[0][2];
    					Tau21k[kind] = U[1][0] + K[1][0];
    					Tau22k[kind] = U[1][1] + K[1][1];
    					Tau23k[kind] = U[1][2] + K[1][2];
    					Tau31k[kind] = U[2][0] + K[2][0];
    					Tau32k[kind] = U[2][1] + K[2][1];
    					Tau33k[kind] = U[2][2] + K[2][2];
                    }
                    else
#endif
                    {
                        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
    						T[i][j] = fftcomplex_t(0.0, 0.0);
    						for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
    							T[i][j] += spectral->C[i][j][k][l]*(U[k][l] + K[k][l]);
    						}
    					}

    					Tau11k[kind] = T[0][0];
    					Tau12k[kind] = T[0][1];
    					Tau13k[kind] = T[0][2];
    					Tau22k[kind] = T[1][1];
    					Tau23k[kind] = T[1][2];
    					Tau33k[kind] = T[2][2];
                    }
				}
			}
		}

		FFT3DBackward((fftw_complex*)Tau11k, (fftw_complex*)Tau11k, N, N, N);
		FFT3DBackward((fftw_complex*)Tau12k, (fftw_complex*)Tau12k, N, N, N);
		FFT3DBackward((fftw_complex*)Tau13k, (fftw_complex*)Tau13k, N, N, N);
		FFT3DBackward((fftw_complex*)Tau22k, (fftw_complex*)Tau22k, N, N, N);
		FFT3DBackward((fftw_complex*)Tau23k, (fftw_complex*)Tau23k, N, N, N);
		FFT3DBackward((fftw_complex*)Tau33k, (fftw_complex*)Tau33k, N, N, N);
#if DISP_GRAD
        if (param->localRotation) {
            FFT3DBackward((fftw_complex*)Tau21k, (fftw_complex*)Tau21k, N, N, N);
            FFT3DBackward((fftw_complex*)Tau31k, (fftw_complex*)Tau31k, N, N, N);
            FFT3DBackward((fftw_complex*)Tau32k, (fftw_complex*)Tau32k, N, N, N);
        }
#endif

		int N3 = N*N*N;
		for (int kx = 0; kx < local_Nx; kx++) {
			for (int ky = 0; ky < N; ky++) {
				for (int kz = 0; kz < N; kz++) {
					int kind = kz+ky*N+kx*N*N;
                    //int kind_global = k+j*N+(i+local_kx_start)*N*N;
                    int kind_local = kz+N*ky+N*N*((myDomain-Dmin)*local_Nx+kx);
#if DISP_GRAD
                    if (param->localRotation) {

                        double G[3][3];
                        G[0][0] = (Tau11k[kind].re)/N3;
    					G[0][1] = (Tau12k[kind].re)/N3;
    					G[0][2] = (Tau13k[kind].re)/N3;
    					G[1][0] = (Tau21k[kind].re)/N3;
    					G[1][1] = (Tau22k[kind].re)/N3;
    					G[1][2] = (Tau23k[kind].re)/N3;
    					G[2][0] = (Tau31k[kind].re)/N3;
    					G[2][1] = (Tau32k[kind].re)/N3;
    					G[2][2] = (Tau33k[kind].re)/N3;

                        double S[3][3];
                        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
    						S[i][j] = 0.0;
    						for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
    							S[i][j] += spectral->C[i][j][k][l]*G[k][l];
    						}
    					}

                        spectral->stress[0][kind_local] = S[0][0];
    					spectral->stress[1][kind_local] = S[1][1];
    					spectral->stress[2][kind_local] = S[2][2];
    					spectral->stress[5][kind_local] = S[0][1];
    					spectral->stress[4][kind_local] = S[0][2];
    					spectral->stress[3][kind_local] = S[1][2];

                        spectral->dispgrad[0][kind_local] = G[0][0];
    					spectral->dispgrad[1][kind_local] = G[0][1];
    					spectral->dispgrad[2][kind_local] = G[0][2];
    					spectral->dispgrad[3][kind_local] = G[1][0];
    					spectral->dispgrad[4][kind_local] = G[1][1];
    					spectral->dispgrad[5][kind_local] = G[1][2];
                        spectral->dispgrad[6][kind_local] = G[2][0];
    					spectral->dispgrad[7][kind_local] = G[2][1];
    					spectral->dispgrad[8][kind_local] = G[2][2];
                    }
                    else
#endif
                    {
                        spectral->stress[0][kind_local] = (Tau11k[kind].re)/N3;
    					spectral->stress[1][kind_local] = (Tau22k[kind].re)/N3;
    					spectral->stress[2][kind_local] = (Tau33k[kind].re)/N3;
    					spectral->stress[5][kind_local] = (Tau12k[kind].re)/N3;
    					spectral->stress[4][kind_local] = (Tau13k[kind].re)/N3;
    					spectral->stress[3][kind_local] = (Tau23k[kind].re)/N3;
                    }

				}
			}
		}

#ifdef PARALLEL

    #if 0
        /* Check everything works - OK */
        MPI_Status  status;

        if (myDomain == 0) {

            double *s0 = (double*)malloc(N*N*N*sizeof(double));
            double *s1 = (double*)malloc(N*N*N*sizeof(double));
            double *s2 = (double*)malloc(N*N*N*sizeof(double));
            double *s3 = (double*)malloc(N*N*N*sizeof(double));
            double *s4 = (double*)malloc(N*N*N*sizeof(double));
            double *s5 = (double*)malloc(N*N*N*sizeof(double));

            memcpy(s0, &spectral->stress[0][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx*sizeof(double));
            memcpy(s1, &spectral->stress[1][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx*sizeof(double));
            memcpy(s2, &spectral->stress[2][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx*sizeof(double));
            memcpy(s3, &spectral->stress[3][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx*sizeof(double));
            memcpy(s4, &spectral->stress[4][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx*sizeof(double));
            memcpy(s5, &spectral->stress[5][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx*sizeof(double));

            for (i = 1; i < numDomains; i++) {
                MPI_Recv(&s0[i*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&s1[i*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&s2[i*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(&s3[i*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
                MPI_Recv(&s4[i*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status);
                MPI_Recv(&s5[i*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status);
            }

            free(spectral->stress[0]);
            free(spectral->stress[1]);
            free(spectral->stress[2]);
            free(spectral->stress[3]);
            free(spectral->stress[4]);
            free(spectral->stress[5]);

            spectral->stress[0] = s0;
            spectral->stress[1] = s1;
            spectral->stress[2] = s2;
            spectral->stress[3] = s3;
            spectral->stress[4] = s4;
            spectral->stress[5] = s5;

        } else {
            MPI_Send(&spectral->stress[0][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&spectral->stress[1][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&spectral->stress[2][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
            MPI_Send(&spectral->stress[3][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            MPI_Send(&spectral->stress[4][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
            MPI_Send(&spectral->stress[5][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
        }
    #else
        MPI_Request request;
        MPI_Status  status;

        for (int d = 0; d < Dspan; d++) {
            if (d == Dspan/2) continue;

            int send = (myDomain + Dspan/2 - d) % numDomains;
            if (send < 0) send += numDomains;

            int recv = (myDomain - Dspan/2 + d) % numDomains;
            if (recv < 0) recv += numDomains;

            for (int k = 0; k < 6; k++) {
                /* Send to MPI task send and receive from MPI task recv */
                MPI_Irecv(&spectral->stress[k][d*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, recv, k, MPI_COMM_WORLD, &request);
                MPI_Send(&spectral->stress[k][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, send, k, MPI_COMM_WORLD);
                MPI_Wait(&request, &status);
            }

#if DISP_GRAD
            if (param->localRotation) {
                for (int k = 0; k < 9; k++) {
                    /* Send to MPI task send and receive from MPI task recv */
                    MPI_Irecv(&spectral->dispgrad[k][d*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, recv, k, MPI_COMM_WORLD, &request);
                    MPI_Send(&spectral->dispgrad[k][(myDomain-Dmin)*N*N*local_Nx], N*N*local_Nx, MPI_DOUBLE, send, k, MPI_COMM_WORLD);
                    MPI_Wait(&request, &status);
                }
            }
#endif
        }
    #endif

#endif

#if 0
        FILE *fp;
        char filename[128];

        sprintf(filename, "stress_%d.dat", myDomain);
        fp = fopen(filename, "w");

        for (kx = 0; kx < local_Nx; kx++)
            for (ky = 0; ky < N; ky++)
                for (kz = 0; kz < N; kz++) {
                    for (k = 0; k < 6; k++)
                        fprintf(fp, "%e ", spectral->stress[k][kz+N*ky+N*N*((myDomain-Dmin)*local_Nx+kx)]);
                    fprintf(fp, "\n");
                }

        fclose(fp);
#endif
		return;
}

/***************************************************************************
 *
 *      Function:    InitCell2Values
 *
 *      Description: A number of values are needed for determining the 'cell2'
 *                   with which a node is to be associated during a cycle.
 *                   Some of these values remain the same throughout the
 *                   entire simulation.  Rather than recalculate these
 *                   'static' values each cycle, we'll calculate them once
 *                   in this function which is only invoked the first time
 *                   nodes must be sorted onto their 'cell2' queues.
 *
 ***************************************************************************/
static void InitCell2Values(Home_t *home,
                            int *cell2PerCellX,
                            int *cell2PerCellY,
                            int *cell2PerCellZ,
                            real8 *cell2SizeX,
                            real8 *cell2SizeY,
                            real8 *cell2SizeZ,
                            real8 *cell2MinX,
                            real8 *cell2MinY,
                            real8 *cell2MinZ,
                            int *numBodyCell2X,
                            int *numBodyCell2Y,
                            int *numBodyCell2Z)
{
        real8   maxSeparation;
        real8   Lx, Ly, Lz;
        real8   probMinX, probMinY, probMinZ;
        real8   probMaxX, probMaxY, probMaxZ;
        real8   cellSizeX, cellSizeY, cellSizeZ;
        Param_t *param;
        Spectral_t *spectral;

        param = home->param;
        spectral = home->spectral;

/*
 *      Get the problem dimensions
 */
        probMinX = param->minSideX; probMaxX = param->maxSideX;
        probMinY = param->minSideY; probMaxY = param->maxSideY;
        probMinZ = param->minSideZ; probMaxZ = param->maxSideZ;

        Lx = probMaxX - probMinX;
        Ly = probMaxY - probMinY;
        Lz = probMaxZ - probMinZ;

        cellSizeX = Lx / param->nXcells;
        cellSizeY = Ly / param->nYcells;
        cellSizeZ = Lz / param->nZcells;

        *cell2MinX = probMinX - cellSizeX;
        *cell2MinY = probMinY - cellSizeY;
        *cell2MinZ = probMinZ - cellSizeZ;

/*
 *      Maximum separation for shrt-range neighboring segments
 */
        maxSeparation = spectral->rcnei;

/*
 *      Compute the size of cell2's, such that an integral number fall within
 *      each cell, and the total number of cell2's in the minimum image problem
 *      space
 */
        *cell2PerCellX = (int) floor(cellSizeX/maxSeparation);
        *cell2PerCellY = (int) floor(cellSizeY/maxSeparation);
        *cell2PerCellZ = (int) floor(cellSizeZ/maxSeparation);

        *cell2PerCellX = MIN(*cell2PerCellX, MAXCELL2PERCELL);
        *cell2PerCellY = MIN(*cell2PerCellY, MAXCELL2PERCELL);
        *cell2PerCellZ = MIN(*cell2PerCellZ, MAXCELL2PERCELL);

        *cell2PerCellX = MAX(*cell2PerCellX, 1);
        *cell2PerCellY = MAX(*cell2PerCellY, 1);
        *cell2PerCellZ = MAX(*cell2PerCellZ, 1);

        *cell2SizeX = cellSizeX / *cell2PerCellX;
        *cell2SizeY = cellSizeY / *cell2PerCellY;
        *cell2SizeZ = cellSizeZ / *cell2PerCellZ;

        *numBodyCell2X = *cell2PerCellX * param->nXcells;
        *numBodyCell2Y = *cell2PerCellY * param->nYcells;
        *numBodyCell2Z = *cell2PerCellZ * param->nZcells;

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FindCell2Index
 *
 *-------------------------------------------------------------------------*/
void SortNodesSpectral(Home_t *home)
{
        int            i, j, k;
        int            minDomCellX, minDomCellY, minDomCellZ;
        int            maxDomCellX, maxDomCellY, maxDomCellZ;
        int            minDomCell2X, minDomCell2Y, minDomCell2Z;
        int            maxDomCell2X, maxDomCell2Y, maxDomCell2Z;
        int            maxNativeCellX, maxNativeCellY, maxNativeCellZ;
        int            minNativeCellX, minNativeCellY, minNativeCellZ;
        int            maxNativeCell2X, maxNativeCell2Y, maxNativeCell2Z;
        int            minNativeCell2X, minNativeCell2Y, minNativeCell2Z;
        int            numCell2X, numCell2Y, numCell2Z, numCell2s;
        int            nextQent = 0;
        Param_t        *param;
        static int     initDone = 0;
        static int     cell2PerCellX, cell2PerCellY, cell2PerCellZ;
        static int     numBodyCell2X, numBodyCell2Y, numBodyCell2Z;
        static real8   cell2SizeX, cell2SizeY, cell2SizeZ;
        static real8   cell2MinX, cell2MinY, cell2MinZ;
        static real8   centerX, centerY, centerZ;


        param = home->param;

/*
 *      The first time into this routine, set some values that won't change
 *      for the remainder of the simulation
 */
        if (!initDone++) {
            InitCell2Values(home,
                            &cell2PerCellX, &cell2PerCellY, &cell2PerCellZ,
                            &cell2SizeX, &cell2SizeY, &cell2SizeZ,
                            &cell2MinX, &cell2MinY, &cell2MinZ,
                            &numBodyCell2X, &numBodyCell2Y, &numBodyCell2Z);
        }

/*
 *      Find the center of the current domain
 */
        centerX = 0.5 * (home->domXmin + home->domXmax);
        centerY = 0.5 * (home->domYmin + home->domYmax);
        centerZ = 0.5 * (home->domZmin + home->domZmax);

/*
 *      Find the current min and max cells for this domain,
 *      including ghost cells
 */
        maxDomCellX = 0;
        maxDomCellY = 0;
        maxDomCellZ = 0;

        minDomCellX = param->nXcells + 2;
        minDomCellY = param->nYcells + 2;
        minDomCellZ = param->nZcells + 2;

/*
 *      NOTE: we *could* thread this loop, but the number of iterations
 *            is likely to be small in most cases AND the overhead for
 *            synchronizing inter-thread access to the values being set
 *            would likely cost more than we would gain from threading
 *            the loop.
 */
        for (i = 0; i < home->nativeCellCount; i++) {
            int     ix, iy, iz;
            int     cellIndex;
            Cell_t  *cell;

            cellIndex = home->cellList[i];
            cell = LookupCell(home, cellIndex);

            ix = cell->xIndex;
            iy = cell->yIndex;
            iz = cell->zIndex;

            if (minDomCellX > ix) minDomCellX = ix;
            if (maxDomCellX < ix) maxDomCellX = ix;
            if (minDomCellY > iy) minDomCellY = iy;
            if (maxDomCellY < iy) maxDomCellY = iy;
            if (minDomCellZ > iz) minDomCellZ = iz;
            if (maxDomCellZ < iz) maxDomCellZ = iz;
        }

        minNativeCellX = minDomCellX;
        minNativeCellY = minDomCellY;
        minNativeCellZ = minDomCellZ;

        maxNativeCellX = maxDomCellX;
        maxNativeCellY = maxDomCellY;
        maxNativeCellZ = maxDomCellZ;

        minDomCellX--;
        minDomCellY--;
        minDomCellZ--;

        maxDomCellX++;
        maxDomCellY++;
        maxDomCellZ++;

/*
 *      Determine the min and max cell2's for this domain
 */
        minDomCell2X = minDomCellX * cell2PerCellX;
        maxDomCell2X = (maxDomCellX + 1) * cell2PerCellX - 1;
        minDomCell2Y = minDomCellY * cell2PerCellY;
        maxDomCell2Y = (maxDomCellY + 1) * cell2PerCellY - 1;
        minDomCell2Z = minDomCellZ * cell2PerCellZ;
        maxDomCell2Z = (maxDomCellZ + 1) * cell2PerCellZ - 1;

        minNativeCell2X = minNativeCellX * cell2PerCellX;
        maxNativeCell2X = (maxNativeCellX+1) * cell2PerCellX - 1;
        minNativeCell2Y = minNativeCellY * cell2PerCellY;
        maxNativeCell2Y = (maxNativeCellY+1) * cell2PerCellY - 1;
        minNativeCell2Z = minNativeCellZ * cell2PerCellZ;
        maxNativeCell2Z = (maxNativeCellZ+1) * cell2PerCellZ - 1;

        numCell2X = maxDomCell2X - minDomCell2X + 1;
        numCell2Y = maxDomCell2Y - minDomCell2Y + 1;
        numCell2Z = maxDomCell2Z - minDomCell2Z + 1;

        numCell2s = numCell2X * numCell2Y * numCell2Z;

/*
 *      Need to save the number cell2 geometry for collision handling
 *      which may need to convert cell2 IDs into the corresponding
 *      X, Y and Z indices of the cell2 in the cell2 grid.
 */
        home->cell2nx = numCell2X;
        home->cell2ny = numCell2Y;
        home->cell2nz = numCell2Z;

/*
 *      Allocate and initialize the array of cell2 queue heads for this domain
 *      and the single array containing the queues of node pointer/array index
 *      pairs.
 */
        home->cell2 = (int *)realloc(home->cell2, numCell2s * sizeof(int));

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < numCell2s; i++) {
            home->cell2[i] = -1;
        }

        home->cell2QentArray = (C2Qent_t *)realloc(home->cell2QentArray,
                                                   home->newNodeKeyPtr *
                                                   sizeof(C2Qent_t));
/*
 *      Loop through all nodes. Queue each node onto one of the cell2's. A node
 *      may fall into more than one cell2 if the domain contains the whole
 *      problem in one or more directions. In this case, queue onto the lowest
 *      numbered image
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            int    baseX, baseY, baseZ;
            int    c2BaseIndex, c2Index;
            int    cell2X, cell2Y, cell2Z;
            int    cell2MinXus, cell2MinYus, cell2MinZus;
            int    cell2Xplus, cell2Yplus, cell2Zplus;
            real8  x, y, z;
            Node_t *node;

            home->cell2QentArray[i].node = (Node_t *)NULL;
            home->cell2QentArray[i].next = -1;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            x = node->x;
            y = node->y;
            z = node->z;

/*
 *          Get the image of the point closest to the center of the
 *          current domain.
 */
            PBCPOSITION(param, centerX, centerY, centerZ, &x, &y, &z);

            cell2X = (int) floor((x - cell2MinX) / cell2SizeX);
            cell2Y = (int) floor((y - cell2MinY) / cell2SizeY);
            cell2Z = (int) floor((z - cell2MinZ) / cell2SizeZ);

/*
 *          If a native node falls outside native cell2 area
 *          force it into the nearest native cell2.
 */
            if (cell2X < minNativeCell2X) {
                cell2X = minNativeCell2X;
            }
            if (cell2X > maxNativeCell2X) {
                cell2X = maxNativeCell2X;
            }

            if (cell2Y < minNativeCell2Y) {
                cell2Y = minNativeCell2Y;
            }
            if (cell2Y > maxNativeCell2Y) {
                cell2Y = maxNativeCell2Y;
            }

            if (cell2Z < minNativeCell2Z) {
                cell2Z = minNativeCell2Z;
            }
            if (cell2Z > maxNativeCell2Z) {
                cell2Z = maxNativeCell2Z;
            }

/*
 *          Find minimum cell2 containing node in the domain range.
 *          This minimum cell will be adjusted to the domain range, rather
 *          than problem range, of cell2 indices.
 *
 *          X direction
 */
            cell2MinXus = cell2X - numBodyCell2X;
            cell2Xplus  = cell2X + numBodyCell2X;

            if (cell2MinXus >= minDomCell2X)
                baseX = cell2MinXus - minDomCell2X;
            else if (cell2X >= minDomCell2X && cell2X <= maxDomCell2X)
                baseX = cell2X - minDomCell2X;
            else if (cell2Xplus <= maxDomCell2X)
                baseX = cell2Xplus - minDomCell2X;
            else
                continue;  /* node moved outside of ghost cells; ignore */

/*
 *          Y direction
 */
            cell2MinYus = cell2Y - numBodyCell2Y;
            cell2Yplus  = cell2Y + numBodyCell2Y;

            if (cell2MinYus >= minDomCell2Y)
                baseY = cell2MinYus - minDomCell2Y;
            else if (cell2Y >= minDomCell2Y && cell2Y <= maxDomCell2Y)
                baseY = cell2Y - minDomCell2Y;
            else if (cell2Yplus <= maxDomCell2Y)
                baseY = cell2Yplus - minDomCell2Y;
            else
                continue; /* node moved outside of ghost cells; ignore */

/*
 *          Z direction
 */
            cell2MinZus = cell2Z - numBodyCell2Z;
            cell2Zplus  = cell2Z + numBodyCell2Z;

            if (cell2MinZus >= minDomCell2Z)
                baseZ = cell2MinZus - minDomCell2Z;
            else if (cell2Z >= minDomCell2Z && cell2Z <= maxDomCell2Z)
                baseZ = cell2Z - minDomCell2Z;
            else if (cell2Zplus <= maxDomCell2Z)
                baseZ = cell2Zplus - minDomCell2Z;
            else
                continue; /* node moved outside of ghost cells; ignore */

/*
 *          Make cell2 index relative to the domain
 */
            cell2X -= minDomCell2X;
            cell2Y -= minDomCell2Y;
            cell2Z -= minDomCell2Z;

/*
 *          Queue the node on the minimum cell2 that it falls into, but save the
 *          ghost cell2 index in the node, since this is where it starts its
 *          neighbor search from.
 *
 *          Use a 'critical' section to prevent multiple threads from updating
 *          the cell2 stuff simultaneously
 */
            c2BaseIndex = EncodeCell2Idx(home, baseX, baseY, baseZ);
            c2Index = EncodeCell2Idx(home, cell2X, cell2Y, cell2Z);

            node->cell2Idx = c2Index;

#ifdef _OPENMP
#pragma omp critical (CRIT_SORT_FOR_COLLISION)
#endif
            {
                node->cell2QentIdx = nextQent;
                home->cell2QentArray[nextQent].node = node;
                home->cell2QentArray[nextQent].next = home->cell2[c2BaseIndex];
                home->cell2[c2BaseIndex] = nextQent;

                nextQent++;

            }  /* end "omp critical" section */
        }

/*
 *      Loop through the cell2s in the domain. If a cell2 is a higher image
 *      of a base cell2 also contained in the domain space, set the higher
 *      image's queue to point to the base image cell2.
 */
        for (i = 0; i < numCell2X; i++) {
            int xIndex, yIndex, zIndex, cell2Index, baseIndex;

            if (i >= numBodyCell2X) {
                xIndex = i - numBodyCell2X;
            } else {
                xIndex = i;
            }

            for (j = 0; j < numCell2Y; j++) {

                if (j >= numBodyCell2Y) {
                    yIndex = j - numBodyCell2Y;
                } else {
                    yIndex = j;
                }

                for (k = 0; k < numCell2Z; k++) {

                    if (k >= numBodyCell2Z) {
                        zIndex = k - numBodyCell2Z;
                    } else {
                        zIndex = k;
                    }

                    baseIndex = EncodeCell2Idx(home, xIndex, yIndex, zIndex);
                    cell2Index = EncodeCell2Idx(home, i, j, k);
                    home->cell2[cell2Index] = home->cell2[baseIndex];
                }
            }
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FindCell2Index
 *
 *-------------------------------------------------------------------------*/
void FindCell2Index(Home_t *home, double px, double py, double pz, int *c2Index)
{
        int            i, j, k;
        int            minDomCellX, minDomCellY, minDomCellZ;
        int            maxDomCellX, maxDomCellY, maxDomCellZ;
        int            minDomCell2X, minDomCell2Y, minDomCell2Z;
        int            maxDomCell2X, maxDomCell2Y, maxDomCell2Z;
        int            maxNativeCellX, maxNativeCellY, maxNativeCellZ;
        int            minNativeCellX, minNativeCellY, minNativeCellZ;
        int            maxNativeCell2X, maxNativeCell2Y, maxNativeCell2Z;
        int            minNativeCell2X, minNativeCell2Y, minNativeCell2Z;
        int            numCell2X, numCell2Y, numCell2Z, numCell2s;
        int            nextQent = 0;
        Param_t        *param;
        static int     initDone = 0;
        static int     cell2PerCellX, cell2PerCellY, cell2PerCellZ;
        static int     numBodyCell2X, numBodyCell2Y, numBodyCell2Z;
        static real8   cell2SizeX, cell2SizeY, cell2SizeZ;
        static real8   cell2MinX, cell2MinY, cell2MinZ;
        static real8   centerX, centerY, centerZ;
        double x, y ,z;
        int    cell2X, cell2Y, cell2Z;
        int    cell2MinXus, cell2MinYus, cell2MinZus;
        int    cell2Xplus, cell2Yplus, cell2Zplus;

        param = home->param;

		x = px;
		y = py;
		z = pz;

/*
 *      The first time into this routine, set some values that won't change
 *      for the remainder of the simulation
 */
        if (!initDone++) {
            InitCell2Values(home,
                            &cell2PerCellX, &cell2PerCellY, &cell2PerCellZ,
                            &cell2SizeX, &cell2SizeY, &cell2SizeZ,
                            &cell2MinX, &cell2MinY, &cell2MinZ,
                            &numBodyCell2X, &numBodyCell2Y, &numBodyCell2Z);
        }

/*
 *      Find the center of the current domain
 */
        centerX = 0.5 * (home->domXmin + home->domXmax);
        centerY = 0.5 * (home->domYmin + home->domYmax);
        centerZ = 0.5 * (home->domZmin + home->domZmax);

/*
 *      Find the current min and max cells for this domain,
 *      including ghost cells
 */
        maxDomCellX = 0;
        maxDomCellY = 0;
        maxDomCellZ = 0;

        minDomCellX = param->nXcells + 2;
        minDomCellY = param->nYcells + 2;
        minDomCellZ = param->nZcells + 2;

/*
 *      NOTE: we *could* thread this loop, but the number of iterations
 *            is likely to be small in most cases AND the overhead for
 *            synchronizing inter-thread access to the values being set
 *            would likely cost more than we would gain from threading
 *            the loop.
 */
        for (i = 0; i < home->nativeCellCount; i++) {
            int     ix, iy, iz;
            int     cellIndex;
            Cell_t  *cell;

            cellIndex = home->cellList[i];
            cell = LookupCell(home, cellIndex);

            ix = cell->xIndex;
            iy = cell->yIndex;
            iz = cell->zIndex;

            if (minDomCellX > ix) minDomCellX = ix;
            if (maxDomCellX < ix) maxDomCellX = ix;
            if (minDomCellY > iy) minDomCellY = iy;
            if (maxDomCellY < iy) maxDomCellY = iy;
            if (minDomCellZ > iz) minDomCellZ = iz;
            if (maxDomCellZ < iz) maxDomCellZ = iz;
        }

        minNativeCellX = minDomCellX;
        minNativeCellY = minDomCellY;
        minNativeCellZ = minDomCellZ;

        maxNativeCellX = maxDomCellX;
        maxNativeCellY = maxDomCellY;
        maxNativeCellZ = maxDomCellZ;

        minDomCellX--;
        minDomCellY--;
        minDomCellZ--;

        maxDomCellX++;
        maxDomCellY++;
        maxDomCellZ++;

/*
 *      Determine the min and max cell2's for this domain
 */
        minDomCell2X = minDomCellX * cell2PerCellX;
        maxDomCell2X = (maxDomCellX + 1) * cell2PerCellX - 1;
        minDomCell2Y = minDomCellY * cell2PerCellY;
        maxDomCell2Y = (maxDomCellY + 1) * cell2PerCellY - 1;
        minDomCell2Z = minDomCellZ * cell2PerCellZ;
        maxDomCell2Z = (maxDomCellZ + 1) * cell2PerCellZ - 1;

        minNativeCell2X = minNativeCellX * cell2PerCellX;
        maxNativeCell2X = (maxNativeCellX+1) * cell2PerCellX - 1;
        minNativeCell2Y = minNativeCellY * cell2PerCellY;
        maxNativeCell2Y = (maxNativeCellY+1) * cell2PerCellY - 1;
        minNativeCell2Z = minNativeCellZ * cell2PerCellZ;
        maxNativeCell2Z = (maxNativeCellZ+1) * cell2PerCellZ - 1;

        numCell2X = maxDomCell2X - minDomCell2X + 1;
        numCell2Y = maxDomCell2Y - minDomCell2Y + 1;
        numCell2Z = maxDomCell2Z - minDomCell2Z + 1;

        numCell2s = numCell2X * numCell2Y * numCell2Z;

/*
 *      Get the image of the point closest to the center of the
 *      current domain.
 */
        PBCPOSITION(param, centerX, centerY, centerZ, &x, &y, &z);

        cell2X = (int) floor((x - cell2MinX) / cell2SizeX);
        cell2Y = (int) floor((y - cell2MinY) / cell2SizeY);
        cell2Z = (int) floor((z - cell2MinZ) / cell2SizeZ);

/*
 *      If a native node falls outside native cell2 area
 *      force it into the nearest native cell2.
 */
        if (cell2X < minNativeCell2X) {
            cell2X = minNativeCell2X;
        }
        if (cell2X > maxNativeCell2X) {
            cell2X = maxNativeCell2X;
        }

        if (cell2Y < minNativeCell2Y) {
            cell2Y = minNativeCell2Y;
        }
        if (cell2Y > maxNativeCell2Y) {
            cell2Y = maxNativeCell2Y;
        }

        if (cell2Z < minNativeCell2Z) {
            cell2Z = minNativeCell2Z;
        }
        if (cell2Z > maxNativeCell2Z) {
            cell2Z = maxNativeCell2Z;
        }

/*
 *      Make cell2 index relative to the domain
 */
        cell2X -= minDomCell2X;
        cell2Y -= minDomCell2Y;
        cell2Z -= minDomCell2Z;

        *c2Index = EncodeCell2Idx(home, cell2X, cell2Y, cell2Z);

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FindNeighborSegments
 *
 *-------------------------------------------------------------------------*/
void FindNeighborSegments(Home_t *home, double px, double py, double pz,
                          int *segsListCnt, int *segsListSize, NeighborSeg_t **segsList)
{
		int     i, j, cell2Index, cx, cy, cz;
		int     cell2X, cell2Y, cell2Z;
		int     nbrCell2Index, nextIndex, arm12;
		double  x1, y1, z1, x2, y2, z2;
		double  vec[3];
		Node_t  *node1, *node2;
		Param_t *param;

		param = home->param;

		*segsListCnt = 0;

#if 1
/*
 * 		Return all segments in the volume
 */
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;

			x1 = node1->x;
			y1 = node1->y;
			z1 = node1->z;

			PBCPOSITION(param, px, py, pz, &x1, &y1, &z1);

			for (j = 0; j < node1->numNbrs; j++) {

				node2 = GetNeighborNode(home, node1, j);
				if (OrderNodes(node1, node2) > 0) continue;

				x2 = node2->x;
				y2 = node2->y;
				z2 = node2->z;

				PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);

				if (*segsListCnt == *segsListSize) {
					*segsListSize += 1000;
					*segsList = (NeighborSeg_t *)realloc(*segsList,
								 sizeof(NeighborSeg_t) * (*segsListSize));
				}

				(*segsList)[*segsListCnt].node1 = node1;
				(*segsList)[*segsListCnt].node2 = node2;
				(*segsList)[*segsListCnt].p1x = x1;
				(*segsList)[*segsListCnt].p1y = y1;
				(*segsList)[*segsListCnt].p1z = z1;
				(*segsList)[*segsListCnt].p2x = x2;
				(*segsList)[*segsListCnt].p2y = y2;
				(*segsList)[*segsListCnt].p2z = z2;
				(*segsList)[*segsListCnt].bx = node1->burgX[j];
				(*segsList)[*segsListCnt].by = node1->burgY[j];
				(*segsList)[*segsListCnt].bz = node1->burgZ[j];

				*segsListCnt += 1;
			}
		}

#else
/*
 * 		Return segments in neighboring cells
 */
		FindCell2Index(home, px, py, pz, &cell2Index);
		if (cell2Index < 0)
			Fatal("cell2Index < 0");

/*
 * 		Loop through all cell2s neighboring the node.  Only
 * 		nodes in these neighboring cell2s are candidates.
 */
		DecodeCell2Idx(home, cell2Index, &cell2X, &cell2Y, &cell2Z);

		for (cx = cell2X - 1; cx <= cell2X + 1; cx++) {
		 for (cy = cell2Y - 1; cy <= cell2Y + 1; cy++) {
		  for (cz = cell2Z - 1; cz <= cell2Z + 1; cz++) {

			nbrCell2Index = EncodeCell2Idx(home, cx, cy, cz);

/*
 *			Loop though all nodes in the neighbor cell2
 */
			nextIndex = home->cell2[nbrCell2Index];

			while (nextIndex >= 0) {

				node1 = home->cell2QentArray[nextIndex].node;
				nextIndex = home->cell2QentArray[nextIndex].next;

				if (node1 == (Node_t *)NULL) continue;

				x1 = node1->x;
				y1 = node1->y;
				z1 = node1->z;

				PBCPOSITION(param, px, py, pz, &x1, &y1, &z1);

/*
 * 				Loop over all arms of node1.
 */
				for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

					node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

					if (node2 == (Node_t *)NULL) continue;

#ifdef PARALLEL
/*
 *					At this point, segment node3/node4 is owned by
 * 					node3.  If node3 is not native to this domain,
 *					the segment may not be used in a collision since
 *					the domain doing to collision must own both
 *					segments.
 */
					if (node3->myTag.domainID != thisDomain) {
						continue;
					}
					if (!DomainOwnsSeg(home, OPCLASS_COLLISION,
								   thisDomain, &node4->myTag)) {
						continue;
					}
#endif

					x2 = node2->x;
					y2 = node2->y;
					z2 = node2->z;

					PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);

/*
 *					It is possible to have a zero-length segment
 *					(created by a previous collision).  If we find
 *					such a segment, do not try to use it in any
 *					subsequent collisions.
 */
					vec[X] = x2 - x1;
					vec[Y] = y2 - y1;
					vec[Z] = z2 - z1;

					if (DotProduct(vec, vec) < 1.0e-20) {
						continue;
					}

/*
 * 					Check if segment if already in the list
 * 					(Need to implement dynamic linked list here)
 */
					int addToList = 1;
					for (i = 0; i < *segsListCnt; i++) {
						if (((*segsList)[i].node1 == node1 && (*segsList)[i].node2 == node2) ||
						    ((*segsList)[i].node1 == node2 && (*segsList)[i].node2 == node1)) {
							addToList = 0;
							break;
						}
					}
/*
 * 					Add to the neighbor list
 */
					if (addToList) {
						if (*segsListCnt == *segsListSize) {
							*segsListSize += 1000;
							*segsList = (NeighborSeg_t *)realloc(*segsList,
										 sizeof(NeighborSeg_t) * (*segsListSize));
						}

						(*segsList)[*segsListCnt].node1 = node1;
						(*segsList)[*segsListCnt].node2 = node2;
						(*segsList)[*segsListCnt].p1x = x1;
						(*segsList)[*segsListCnt].p1y = y1;
						(*segsList)[*segsListCnt].p1z = z1;
						(*segsList)[*segsListCnt].p2x = x2;
						(*segsList)[*segsListCnt].p2y = y2;
						(*segsList)[*segsListCnt].p2z = z2;
						(*segsList)[*segsListCnt].bx = node1->burgX[arm12];
						(*segsList)[*segsListCnt].by = node1->burgY[arm12];
						(*segsList)[*segsListCnt].bz = node1->burgZ[arm12];

						*segsListCnt += 1;
					}

				}  /* Loop over node3 arms */
			}  /* while(nextIndex >= 0) */
		   }  /* Loop over neighboring cell2s */
	      }  /* Loop over neighboring cell2s */
	     }  /* Loop over neighboring cell2s */
#endif

		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SpectralDispGrad
 *
 *-------------------------------------------------------------------------*/
void SpectralDispGrad(Home_t *home, double px, double py, double pz, double G[9])
{
	int      i, j, k, l, N, kind;
	double   p[3], H[3], q[3], xi[3], phi;
	int      g[3];
	double   MU, NU, a;
	Param_t  *param;
	Spectral_t *spectral;

	param = home->param;
	spectral = home->spectral;

    if (!param->localRotation)
        Fatal("SpectralDispGrad function only available for localRotation=1");

#ifdef PARALLEL
    int myDomain = home->myDomain;
    int numDomains = home->numDomains;
    int Dspan = spectral->Dspan;
    int Dmin = myDomain-Dspan/2;
    ptrdiff_t local_Nx = spectral->local_Nx;
#endif

	MU = param->shearModulus;
	NU = param->pois;
	a  = param->rc;

	N = spectral->N;
	H[0] = param->Lx/N;
	H[1] = param->Ly/N;
	H[2] = param->Lz/N;

#ifdef PARALLEL
    double cx = param->minSideX + (myDomain+0.5)*param->Lx/numDomains;
    double sx = ( (param->xBoundType == Periodic) ? (1.0/param->Lx) : 0.0 );

    p[0] = px - rint((px-cx)*sx) * param->Lx;   // (apply PBC correction)
    p[0] -= param->minSideX;
#else
	p[0] = px - param->minSideX;
#endif
    p[1] = py - param->minSideY;
    p[2] = pz - param->minSideZ;

	q[0] = p[0] / H[0] - 0.5;
	q[1] = p[1] / H[1] - 0.5;
	q[2] = p[2] / H[2] - 0.5;

	g[0] = (int)floor(q[0]);
	g[1] = (int)floor(q[1]);
	g[2] = (int)floor(q[2]);

	xi[0] = 2.0*(q[0]-g[0]) - 1.0;
	xi[1] = 2.0*(q[1]-g[1]) - 1.0;
	xi[2] = 2.0*(q[2]-g[2]) - 1.0;

	/* Determine elements for interpolation and apply PBC */
	int ind1d[3][2];
	for (i = 0; i < 2; i++) {
#ifdef PARALLEL
        ind1d[0][i] = (g[0]+i) - Dmin*local_Nx;
        if (ind1d[0][i] < 0 || ind1d[0][i] >= Dspan*local_Nx)
            Fatal("Local disp grad out of Dpsan on proc %d", myDomain);
#else
        ind1d[0][i] = (g[0]+i)%N;
		if (ind1d[0][i] < 0) ind1d[0][i] += N;
#endif
		ind1d[1][i] = (g[1]+i)%N;
		if (ind1d[1][i] < 0) ind1d[1][i] += N;
		ind1d[2][i] = (g[2]+i)%N;
		if (ind1d[2][i] < 0) ind1d[2][i] += N;
	}

	// 1d shape functions
	double phi1d[3][2];
	for (i = 0; i < 3; i++) {
		phi1d[i][0] = 0.5*(1.0-xi[i]);
		phi1d[i][1] = 0.5*(1.0+xi[i]);
	}

	// 3d shape functions and indices
	for (l = 0; l < 9; l++) {
		G[l] = 0.0;
		int a=0;
		for (k = 0; k < 2; k++) {
			for (j = 0; j < 2; j++) {
				for (i = 0; i < 2; i++) {
					phi = phi1d[0][i]*phi1d[1][j]*phi1d[2][k];
                    kind = ind1d[2][k]+ind1d[1][j]*N+ind1d[0][i]*N*N;
					G[l] += phi * spectral->dispgrad[l][kind];
				}
			}
		}
	}
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SpectralStress
 *
 *-------------------------------------------------------------------------*/
void SpectralStress(Home_t *home, int includeNei,
                    double px, double py, double pz, double Stotal[6])
{
	int      i, j, k, l, N, kind;
	double   p[3], H[3], q[3], xi[3], phi;
	int      g[3];
	double   MU, NU, a;
	Param_t  *param;
	Spectral_t *spectral;

	param = home->param;
	spectral = home->spectral;

#ifdef PARALLEL
    int myDomain = home->myDomain;
    int numDomains = home->numDomains;
    int Dspan = spectral->Dspan;
    int Dmin = myDomain-Dspan/2;
    ptrdiff_t local_Nx = spectral->local_Nx;
#endif

	MU = param->shearModulus;
	NU = param->pois;
	a  = param->rc;

	N = spectral->N;
	H[0] = param->Lx/N;
	H[1] = param->Ly/N;
	H[2] = param->Lz/N;

#ifdef PARALLEL
    double cx = param->minSideX + (myDomain+0.5)*param->Lx/numDomains;
    double sx = ( (param->xBoundType == Periodic) ? 1.0/param->Lx : 0.0 );
    p[0] = px - rint((px-cx)*sx) * param->Lx;   // apply PBC correction <x>
    p[0] -= param->minSideX;
#else
	p[0] = px - param->minSideX;
#endif
    p[1] = py - param->minSideY;
    p[2] = pz - param->minSideZ;

	q[0] = p[0] / H[0] - 0.5;
	q[1] = p[1] / H[1] - 0.5;
	q[2] = p[2] / H[2] - 0.5;

	g[0] = (int)floor(q[0]);
	g[1] = (int)floor(q[1]);
	g[2] = (int)floor(q[2]);

	xi[0] = 2.0*(q[0]-g[0]) - 1.0;
	xi[1] = 2.0*(q[1]-g[1]) - 1.0;
	xi[2] = 2.0*(q[2]-g[2]) - 1.0;

	/* Determine elements for interpolation and apply PBC */
    int skip = 0;
	int ind1d[3][2];
	for (i = 0; i < 2; i++) {
#ifdef PARALLEL
        ind1d[0][i] = (g[0]+i) - Dmin*local_Nx;
        if (ind1d[0][i] < 0 || ind1d[0][i] >= Dspan*local_Nx) {
#if 0
            Fatal("Local stress out of Dpsan on proc %d", myDomain);
#else
            skip = 1;
            printf("Warning: Local stress out of Dpsan on proc %d\n", myDomain);
#endif
        }
#else
        ind1d[0][i] = (g[0]+i)%N;
		if (ind1d[0][i] < 0) ind1d[0][i] += N;
#endif
		ind1d[1][i] = (g[1]+i)%N;
		if (ind1d[1][i] < 0) ind1d[1][i] += N;
		ind1d[2][i] = (g[2]+i)%N;
		if (ind1d[2][i] < 0) ind1d[2][i] += N;
	}

	// 1d shape functions
	double phi1d[3][2];
	for (i = 0; i < 3; i++) {
		phi1d[i][0] = 0.5*(1.0-xi[i]);
		phi1d[i][1] = 0.5*(1.0+xi[i]);
	}

	// 3d shape functions and indices
	for (l = 0; l < 6; l++) {
		Stotal[l] = 0.0;
		if (skip) continue;
		for (k = 0; k < 2; k++) {
			for (j = 0; j < 2; j++) {
				for (i = 0; i < 2; i++) {
					phi = phi1d[0][i]*phi1d[1][j]*phi1d[2][k];
                    kind = ind1d[2][k]+ind1d[1][j]*N+ind1d[0][i]*N*N;
					Stotal[l] += phi * spectral->stress[l][kind];
				}
			}
		}
	}

#ifndef PARALLEL
/*
 * 	Extra stress contribution from segment portions within h/2
 */
	if (includeNei == INCLUDE_NEI && spectral->rcnei > 0.0) {

		int numSegs, segsListSize;
		NeighborSeg_t *segsList;

		segsListSize = 10000;
		segsList = (NeighborSeg_t *)calloc(segsListSize, sizeof(NeighborSeg_t));

		FindNeighborSegments(home, px, py, pz, &numSegs, &segsListSize, &segsList);

		p[0] = px;
		p[1] = py;
		p[2] = pz;

		double r, r2;
		r = spectral->rcnei;
		r2 = r*r;

		//printf("  numNeigbhor = %d, r = %e\n",numSegs,r);

		int numAdd = 0;

		for (i = 0; i < numSegs; i++) {

			int    addStress = 0;
			int    in1, in2;
			double p1[3], p2[3], p1in[3], p2in[3];
			double mu1, mu2, mutmp;

			p1[0] = segsList[i].p1x;
			p1[1] = segsList[i].p1y;
			p1[2] = segsList[i].p1z;

			p2[0] = segsList[i].p2x;
			p2[1] = segsList[i].p2y;
			p2[2] = segsList[i].p2z;

/*
 * 		   Include whole segment if the minimum distance is closer than rcnei
 */
			double d0, d1, d2;
			double dt, t[3], pp1[3], pp2[3], u[3];

			t[0] = p2[0] - p1[0];
			t[1] = p2[1] - p1[1];
			t[2] = p2[2] - p1[2];
			dt = Normal(t);

			pp1[0] = p1[0] - p[0];
			pp1[1] = p1[1] - p[1];
			pp1[2] = p1[2] - p[2];

			pp2[0] = p2[0] - p[0];
			pp2[1] = p2[1] - p[1];
			pp2[2] = p2[2] - p[2];

			CrossVector(pp1, t, u);
			d0 = Normal(u)/dt;
			d1 = Normal(pp1);
			d2 = Normal(pp2);

			d0 = MIN(d0,d1);
			d0 = MIN(d0,d2);

			if (dt > 1.e-10 && d0 < r) {
				VECTOR_COPY(p1in, p1);
				VECTOR_COPY(p2in, p2);
				addStress = 1;
			}

			if (addStress) {
				//printf("  add stress segment %d\n", i);
				numAdd++;

				double bX, bY, bZ;
				bX = segsList[i].bx;
				bY = segsList[i].by;
				bZ = segsList[i].bz;

				double Sseg0[6], Sseg1[6], Sseg11[6], Sseg12[6];

				StressDueToSeg(home, px, py, pz, p1in[0], p1in[1], p1in[2],
							   p2in[0], p2in[1], p2in[2], bX, bY, bZ,
							   spectral->rcgrid, MU, NU, Sseg0);

#if SINGLY_CONVOLUTED
                /* Compute the singly-convoluted analytical stress */
                StressDueToSeg(home, px, py, pz, p1in[0], p1in[1], p1in[2],
                               p2in[0], p2in[1], p2in[2], bX, bY, bZ,
                               0.9038*a, MU, NU, Sseg11);

                StressDueToSeg(home, px, py, pz, p1in[0], p1in[1], p1in[2],
                               p2in[0], p2in[1], p2in[2], bX, bY, bZ,
                               0.5451*a, MU, NU, Sseg12);

                for (l = 0; l < 6; l++) {
                    Stotal[l] += (0.3425*Sseg11[l] + 0.6575*Sseg12[l] - Sseg0[l]);
                }
#else
				StressDueToSeg(home, px, py, pz, p1in[0], p1in[1], p1in[2],
							   p2in[0], p2in[1], p2in[2], bX, bY, bZ,
							   a, MU, NU, Sseg1);

				for (l = 0; l < 6; l++) {
					Stotal[l] += (Sseg1[l] - Sseg0[l]);
				}
#endif
			}
		}

		//printf("  numNei stress = %d\n",numAdd);

		free(segsList);
	}
#else
    if (includeNei == INCLUDE_NEI)
        Fatal("SpectralStress function does not support INCLUDE_NEI in parallel");
#endif

}

/*---------------------------------------------------------------------------
 *
 *      Function:     SpectralSigbForce
 *
 *-------------------------------------------------------------------------*/
void SpectralSigbForce(Home_t *home, double p1[3], double p2[3],
                       double burg[3], double f1[3], double f2[3])
{
		int      i, numPoints;
		real8    *positions, *weights;
		real8    temp, mult1, mult2;
		real8    sigbx, sigby, sigbz;
		real8    fLinvx, fLinvy, fLinvz;
		real8    pmidx, pmidy, pmidz;
		real8    pspanx, pspany, pspanz;
		real8    pos, p[3], S[6];
		Param_t  *param;
		Spectral_t *spectral;

		param  = home->param;
		spectral = home->spectral;

		numPoints = spectral->intNumPoints;
		positions = spectral->glPositions;
		weights   = spectral->glWeights;

		f1[0] = 0.0;
		f1[1] = 0.0;
		f1[2] = 0.0;

		f2[0] = 0.0;
		f2[1] = 0.0;
		f2[2] = 0.0;

		pmidx  = 0.5 * (p2[X]+p1[X]);
		pmidy  = 0.5 * (p2[Y]+p1[Y]);
		pmidz  = 0.5 * (p2[Z]+p1[Z]);

		pspanx = 0.5 * (p2[X]-p1[X]);
		pspany = 0.5 * (p2[Y]-p1[Y]);
		pspanz = 0.5 * (p2[Z]-p1[Z]);


		for (i = 0; i < numPoints; i++) {

			pos = positions[i];
			p[X] = pmidx+pspanx*pos;
			p[Y] = pmidy+pspany*pos;
			p[Z] = pmidz+pspanz*pos;

			SpectralStress(home, NO_INCLUDE_NEI, p[X], p[Y], p[Z], S);

            sigbx = S[0]*burg[X] + S[5]*burg[Y] + S[4]*burg[Z];
			sigby = S[5]*burg[X] + S[1]*burg[Y] + S[3]*burg[Z];
			sigbz = S[4]*burg[X] + S[3]*burg[Y] + S[2]*burg[Z];

			fLinvx = (sigby*pspanz-sigbz*pspany);
			fLinvy = (sigbz*pspanx-sigbx*pspanz);
			fLinvz = (sigbx*pspany-sigby*pspanx);

			temp = weights[i]*positions[i];
			mult1 = weights[i]+temp;

			f2[0] += fLinvx*mult1;
			f2[1] += fLinvy*mult1;
			f2[2] += fLinvz*mult1;

			mult2 = weights[i]-temp;

			f1[0] += fLinvx*mult2;
			f1[1] += fLinvy*mult2;
			f1[2] += fLinvz*mult2;
		}

		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     SpectralForceOneSeg
 *
 *-------------------------------------------------------------------------*/
void SpectralForceOneSeg(Home_t *home, Node_t *node1, Node_t *node2,
                         real8 f1[3], real8 f2[3])
{
        int       armID;
        real8     p1[3], p2[3], burg[3];
        real8     p1f[3], p2f[3];
        Param_t   *param;

        param = home->param;

        armID = GetArmID(node1, node2);

        p1[X] = node1->x;
        p1[Y] = node1->y;
        p1[Z] = node1->z;

        p2[X] = node2->x;
        p2[Y] = node2->y;
        p2[Z] = node2->z;

        burg[X] = node1->burgX[armID];
        burg[Y] = node1->burgY[armID];
        burg[Z] = node1->burgZ[armID];

        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);

        SpectralSigbForce(home, p1, p2, burg, p1f, p2f);

/*
 *      Return the arm-specific forces for both the nodes of the segment,
 *      the total nodal forces will be set to the sum of the arm forces later.
 */
        f1[X] = p1f[X];
        f1[Y] = p1f[Y];
        f1[Z] = p1f[Z];

        f2[X] = p2f[X];
        f2[Y] = p2f[Y];
        f2[Z] = p2f[Z];

        return;
}

#if 0
/*---------------------------------------------------------------------------
 *
 *      Function:       SpectralNodeForce
 *
 *-------------------------------------------------------------------------*/
void SpectralNodeForce(Home_t *home)
{
		int     i, j, armID12, armID21;
		double  p1[3], p2[3], dx, dy, dz;
		double  burg[3], f1[3], f2[3];
		Node_t  *node1, *node2;
		Param_t *param;

		param = home->param;

        ZeroNodeForces(home, FULL);
        SortNodesSpectral(home);
        FFTCellCharge(home);


		for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            p1[X] = node1->x;
            p1[Y] = node1->y;
            p1[Z] = node1->z;

            for (j = 0; j < node1->numNbrs; j++) {

                node2 = GetNeighborNode(home, node1, j);
/*
 *              Insures node1 is the node with the lower tag
 */
                if (OrderNodes(node1, node2) >= 0) {
                    continue;
                }

                dx = node2->x - p1[X];
                dy = node2->y - p1[Y];
                dz = node2->z - p1[Z];

                ZImage(param, &dx, &dy, &dz);

/*
 *             It is possible to have a zero-length segment (created by
 *             collision handling).  If we find such a segment, there will
 *             be no forces on the segment, so just skip to the next segment.
 */
                if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                    continue;
                }

                p2[X] = p1[X] + dx;
                p2[Y] = p1[Y] + dy;
                p2[Z] = p1[Z] + dz;

                burg[X] = node1->burgX[j];
                burg[Y] = node1->burgY[j];
                burg[Z] = node1->burgZ[j];

                armID12 = j;
                armID21 = GetArmID(node2, node1);

/*
 *              Add in force due to long-range spectral stress
 */
				SpectralSigbForce(home, p1, p2, burg, f1, f2);

				AddtoArmForce(node1, armID12, f1);
				AddtoArmForce(node2, armID21, f2);
			}
		}
#if 1 // now done in NodeForceList2
/*
 *      We should now have updated forces for nodes/segments
 *      so now do a quick loop through all local nodes and set
 *      the nodes' total forces to the sum of their arms' forces.
 */
		//TimerStart(home, ADD_FORCE);

		#ifdef _OPENMP
	/*	#pragma omp parallel for default(none) schedule(dynamic,1) \ */
		#pragma omp parallel for default(none) schedule(static) \
			private(i , j , node1) \
			shared (home  )
		#endif
		for (i = 0; i < home->newNodeKeyPtr; i++) {

			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			//if (node1->subgroup == 0) continue;

			node1->fX = 0.0;
			node1->fY = 0.0;
			node1->fZ = 0.0;

			for (j = 0; j < node1->numNbrs; j++) {
				node1->fX += node1->armfx[j];
				node1->fY += node1->armfy[j];
				node1->fZ += node1->armfz[j];
			}
		}
		//TimerStop(home, ADD_FORCE);
#endif

		return;
}
#endif

/*---------------------------------------------------------------------------
 *
 *      Function:     FFTCellCharge
 *
 *-------------------------------------------------------------------------*/
void FFTCellCharge(Home_t *home)
{
        FFTAssignSlice(home);
        AlphaTensor(home);
        StressFromAlpha(home);

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FinishSpectral
 *
 *-------------------------------------------------------------------------*/
void FinishSpectral(Home_t *home)
{
		int i, j, k, l, N;
        Param_t *param;
		Spectral_t *spectral;

        param = home->param;
		spectral = home->spectral;

        if (!param->FFTenabled) return;

        for (i = 0; i < 9; i++)
            free(spectral->fields[i]);

#if NS_DISTRIB
        free(spectral->wk);
#endif

		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     InitSpectralCore
 *
 *-------------------------------------------------------------------------*/
void InitSpectralCore(Home_t *home)
{
        int N;
        Param_t *param;
        Spectral_t *spectral;

        param = home->param;
		spectral = home->spectral;
		N = spectral->N;

        double H[3], Hmax;
        H[0] = param->Lx/N;
        H[1] = param->Ly/N;
        H[2] = param->Lz/N;
        Hmax = MAX(MAX(H[0], H[1]), H[2]);

#if SINGLY_CONVOLUTED
        spectral->rcgrid = 4.0*Hmax;
#else
        spectral->rcgrid = 2.0*Hmax;
#endif

        // If we are using line-tension with FFT, set the core radius
        // to the effective grid core radius and recalculate core energy

        if (!param->elasticinteraction) {
            param->rc = spectral->rcgrid;
            param->Ecore = (param->shearModulus / (4*M_PI)) *
                           log(param->rc/0.1);
            param->TensionFactor = 2.0 * param->Ecore / param->shearModulus;
        }

        if (param->rc >= spectral->rcgrid) {
            spectral->rcgrid = param->rc;
            spectral->rcnei = 0.0;
        } else {
            spectral->rcnei = 4.0*spectral->rcgrid;
        }

#if NS_DISTRIB
        /* Non singular distribution kernel */
        int kmax = N/2 + N % 2;

		if (fabs(H[0]-H[1]) > 1e-10 || fabs(H[0]-H[2]) > 1e-10)
			printf("WARNING: non-cubic volume = %e %e %e\n", param->Lx, param->Ly, param->Lz);

		double aw = spectral->rcgrid;
        double aw2 = aw*aw;
        double aw4 = aw2*aw2;

#ifdef PARALLEL
        ptrdiff_t local_Nx = spectral->local_Nx;
        ptrdiff_t local_kx_start = spectral->local_kx_start;
#else
        int local_Nx = spectral->local_Nx;
        int local_kx_start = spectral->local_kx_start;
#endif

        int kxx, kyy, kzz, kind;
		double wtot = 0.0;
		double x, y, z, r2;
		for (int i = 0; i < local_Nx; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < N; k++) {
                    kxx = i + local_kx_start; kyy = j; kzz = k;
					if (kxx >= kmax) kxx = N - kxx;
					if (kyy >= kmax) kyy = N - kyy;
					if (kzz >= kmax) kzz = N - kzz;
					x = kxx*H[0];
					y = kyy*H[1];
					z = kzz*H[2];
					r2 = x*x + y*y + z*z;
                    kind = k+j*N+i*N*N;
					spectral->wk[kind] = fftcomplex_t(15.0*aw4/8.0/M_PI/pow(r2+aw2, 3.5), 0.0);
					wtot += spectral->wk[kind].re;
				}
			}
		}
#ifdef PARALLEL
        MPI_Allreduce(MPI_IN_PLACE, &wtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

		for (int i = 0; i < local_Nx; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < N; k++) {
					kind = k+j*N+i*N*N;
					spectral->wk[kind].re /= wtot;
				}
			}
		}

        FFT3DForward((fftw_complex*)spectral->wk, (fftw_complex*)spectral->wk, N, N, N);
#endif

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     AllocateSpectral
 *
 *-------------------------------------------------------------------------*/
void AllocateSpectral(Home_t *home)
{
		int i, j, k, l, N, Nx;
        Param_t *param;
		Spectral_t *spectral;

        param = home->param;
		spectral = home->spectral;
		N = spectral->N;

#ifdef PARALLEL
        if (N % home->numDomains != 0)
            Fatal("The number of MPI tasks must be a divisor of the FFT grid size");
        if (param->nXdoms != home->numDomains)
            Fatal("The domain must only be partitionned along the x direction when using the FFT");

        fftw_mpi_init();
        ptrdiff_t local_Nx, local_kx_start;
        fftw_mpi_local_size_3d(N, N, N, MPI_COMM_WORLD, &local_Nx, &local_kx_start);
        spectral->Dspan = 1;
#else
        int local_Nx = N;
        int local_kx_start = 0;
#endif

        spectral->local_Nx = local_Nx;
        spectral->local_kx_start = local_kx_start;

        for (i = 0; i < 9; i++) {
			spectral->fields[i] = (fftcomplex_t*)malloc(local_Nx*N*N*sizeof(fftcomplex_t));
			if (spectral->fields[i] == NULL)
                Fatal("malloc error for spectral fields array");

            spectral->alpha[i] = spectral->fields[i];
		}

        int stressptr[6] = {0, 4, 8, 5, 2, 1};
        for (i = 0; i < 6; i++)
            spectral->stress[i] = (double*)spectral->fields[stressptr[i]];

        if (param->localRotation) {
#if DISP_GRAD
            // Allocate displacement gradient
            for (i = 0; i < 9; i++) {
    			spectral->dispgrad[i] = (double*)malloc(local_Nx*N*N*sizeof(double));
    			if (spectral->dispgrad[i] == NULL)
                    Fatal("malloc error for displacement gradient field array");
    		}
#else
            Fatal("The code must be compiled with DISP_GRAD to enable localRotation");
#endif
        }

#if NS_DISTRIB
        /* Non singular distribution kernel */
		spectral->wk = (fftcomplex_t*)malloc(sizeof(fftcomplex_t)*local_Nx*N*N);
#endif

        spectral->intNumPoints = 16; //128;
		spectral->glPositions = (double *)malloc(spectral->intNumPoints * sizeof(double));
		spectral->glWeights   = (double *)malloc(spectral->intNumPoints * sizeof(double));

		GaussQuadCoeffAll(spectral->intNumPoints, spectral->glPositions, spectral->glWeights);

		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     InitSpectralModulus
 *
 *-------------------------------------------------------------------------*/
void InitSpectralModulus(Home_t *home)
{

    Param_t    *param    = home->param;
    Spectral_t *spectral = home->spectral;

#ifndef ANISOTROPIC
    real8 MU = param->shearModulus;
    real8 NU = param->pois;
    real8 LA = 2.0*MU*NU/(1.0-2.0*NU);

    real8 delta[3][3] = { 1, 0, 0,
                          0, 1, 0,
                          0, 0, 1 };

    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    for (int k=0; k<3; k++)
    for (int l=0; l<3; l++) {
       spectral->C[i][j][k][l] =   LA*(delta[i][j]*delta[k][l]) 
                                 + MU*(delta[i][k]*delta[j][l]+delta[i][l]*delta[j][k]);
    }
#else
    SetElasticConstantMatrix4D(param->elasticConstantMatrix, spectral->C);
#endif

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     InitSpectral
 *
 *-------------------------------------------------------------------------*/
void InitSpectral(Home_t *home)
{
	Param_t  *param = home->param;

    if (!param->FFTenabled) return;

	Spectral_t *spectral = (Spectral_t*)calloc(1, sizeof(Spectral_t));

	home->spectral = spectral;
	spectral->N    = param->FFTgrid;

	AllocateSpectral(home);

    InitSpectralCore(home);

	InitSpectralModulus(home);

    if (home->myDomain == 0) {
        printf("***********************************************************\n");
        printf("**** \n");
        printf("**** DDD spectral module (DDD-FFT)\n");
        printf("**** \n");
        printf("***********************************************************\n");
    	printf("Core radius: a = %e, agrid = %e\n", param->rc, spectral->rcgrid);
    	printf("Grid size: FFTgrid = %d x %d x %d\n", spectral->N, spectral->N, spectral->N);
#ifndef PARALLEL
    	printf("Number of dislocation nodes: %d\n", home->newNodeKeyPtr);
#endif
    	printf("\n");
    }

    // Deactivate load balancing? 
    // param->numDLBCycles = 0;
    // param->DLBfreq = 0;
}

#endif
