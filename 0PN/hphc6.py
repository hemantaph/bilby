import numpy as np

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22

#(* Euler-Gamma constant. Appears in 3PN corrections *) 
gamma = 0.5772156649 

class Fn:
    def __init__(self, iota_, beta_, D_, m1_, m2_, f_, f0_, Fp_, Fc_, et0_, phic_, tc_ ):

        self.iota_ = iota_
        self.beta_ = beta_
        self.D_ = D_*Mpc
        self.m1_ = m1_*Mo
        self.m2_ = m2_*Mo
        self.f0_ = f0_
        self.Fp_ = Fp_
        self.Fc_ = Fc_
        self.et0_ = et0_
        self.phic_ = phic_
        self.tc_ = tc_
        self.f_ = f_
        self.M_ = (m1_+m2_)*Mo
        self.eta_ = (m1_*m2_)/((m1_+m2_)**2)
        self.chi_ = f_/f0_
        self.ff_ = C**3/(G*(m1_+m2_)*Mo*np.pi*6**(3/2))
        self.delta_ = (m1_-m2_)*Mo
    
    #defining unit-step function 
    def unitstep(self,lp,ff,f):
        if lp*ff-2*f>=0:
            return(1)
        else:
            return(0)
        
    #it is plausible to take small sine reullt as 0
    #coz sine(pi) is 0, but np.sin(np.pi) is not 0
    #but u cannot the resolve the beta and iota angle between pi and np.pi 
    def sine_(self, x_):
        result_ = np.sin(float(x_)) 
        if abs(result_)<=4.898587196589413e-16:
            result_ = 0.0
            
        return(result_)
    
    def cosine_(self, x_):
        result_ = np.cos(float(x_)) 
        if abs(result_)<=1.8369701987210297e-16:
            result_ = 0.0
            
        return(result_)
    
    def htilde(self):
        iota = self.iota_
        beta = self.beta_ 
        D = self.D_
        m1 = self.m1_
        m2 = self.m2_
        f0 = self.f0_
        Fp = self.Fp_
        Fc = self.Fc_
        et0 = self.et0_
        phic = self.phic_
        tc = self.tc_
        f = self.f_
        M = self.M_
        eta = self.eta_
        chi = self.chi_
        ff = self.ff_
        delta = self.delta_
        sine = self.sine_
        cosine = self.cosine_
        
        #k depends only on l
        l = np.array([1,2,3,4,5,6,7,8,9,10])
        xk = ((G*M*2*np.pi*f)/( C**3 * l))**(2/3)
        
        k = (xk**2 * (27/2 - 7*eta + et0**6 * (1589015535/(5175296*chi**7) - (110496315*eta)/(184832*chi**7) + 498132560585/(1614692352*chi**(19/3)) + (168333575*eta)/(1478656*chi**(19/3)) - 908911888607/(1816528896*chi**(44/9)) + (4154545993*eta)/(4990464*chi**(44/9)) - 892815840919/(1816528896*chi**(38/9)) - (301299733*eta)/(4990464*chi**(38/9)) + 3110697122471/(14532231168*chi**(25/9)) - (12075357445*eta)/(39923712*chi**(25/9)) + 17781607555/(69866496*chi**(19/9)) - (82799465*eta)/(2495232*chi**(19/9))) + et0**4 * ( - (6830363/(153216*chi**(44/9))) + (474967*eta)/(5472*chi**(44/9)) - 268677653/(3983616*chi**(38/9)) - (90671*eta)/(10944*chi**(38/9)) + 18185905/(284544*chi**(25/9)) - (267133*eta)/(2736*chi**(25/9)) + 34967929/(306432*chi**(19/9)) - (162827*eta)/(10944*chi**(19/9))) + et0**2 * (2833/(336*chi**(25/9)) - (197*eta)/(12*chi**(25/9)) + 10523/(336*chi**(19/9)) - (49*eta)/(12*chi**(19/9)))) + xk * (3 + et0**4 * ( - (2411/(304*chi**(38/9))) + 3323/(304*chi**(19/9))) + et0**6 * (1682685/(46208*chi**(19/3)) - 8011753/(138624*chi**(38/9)) + 1689785/(69312*chi**(19/9))) + (3 * et0**2)/chi**(19/9)) + xk**(5/2) * ( - ((377 * et0**2*np.pi*( - 1 + chi))/(24*chi**(28/9))) - (et0**4 * np.pi*(3635788 - 2258257*chi - 3883073*chi**(19/9) + 2505542*chi**(28/9)))/(43776*chi**(47/9)) - (et0**6 * np.pi*( - 331142311890 + 189286331247*chi + 446685572185*chi**(19/9) - 217621452319*chi**(28/9) - 161105816843*chi**(38/9) + 73897677620*chi**(47/9)))/(578893824*chi**(22/3))) + xk**3 * (1/32*(2160 - 5192*eta + 123*np.pi**2*eta + 224*eta**2) + 1/(1016064*chi**(31/9)) * et0**2 * ( - 1193251 - 22282512*eta + 42700560*eta**2 + 89434977*chi**(2/3) - 185795232*eta*chi**(2/3) + 22703856*eta**2*chi**(2/3) + 285923842*chi**(4/3) - 477003408*eta*chi**(4/3) + 12692862*np.pi**2*eta*chi**(4/3) - 3424512*eta**2*chi**(4/3)) + 1/(84325183488*chi**(50/9)) * et0**4 * ( - 4759063292165 + 30348860970792*eta - 38763744107088*eta**2 - 31968879219858*chi**(2/3) + 58318046249232*eta*chi**(2/3) + 7646159215968*eta**2*chi**(2/3) - 31736268138496*chi**(4/3) + 45889671128952*eta*chi**(4/3) - 1650389381550*np.pi**2*eta*chi**(4/3) + 219411602592*eta**2*chi**(4/3) + 2043730880707*chi**(19/9) - 27141856655592*eta*chi**(19/9) + 35073241576464*eta**2*chi**(19/9) + 56262861824610*chi**(25/9) - 93285988260024*eta*chi**(25/9) + 11206302010176*eta**2*chi**(25/9) + 86461368353906*chi**(31/9) - 144242491555344*eta*chi**(31/9) + 3838232618766*np.pi**2*eta*chi**(31/9) - 1035548457216*eta**2*chi**(31/9)) + 1/(3999037501734912*chi**(23/3)) * et0**6 * (3279894679024104960 - 16631859540896308800*eta + 19202503791507079680*eta**2 + 10402025549836075155*chi**(2/3) - 16414748179030197360*eta*chi**(2/3) - 7473844435281260400*eta**2*chi**(2/3) + 6636901757399139267*chi**(4/3) - 6647849008784363016*eta*chi**(4/3) + 341159644201429872*np.pi**2*eta*chi**(4/3) + 623118469564347120*eta**2*chi**(4/3) - 4220015697065866256*chi**(19/9) + 21692272309533291648*eta*chi**(19/9) - 23329096412641114368*eta**2*chi**(19/9) - 17016310487961012648*chi**(25/9) + 26221595034201700224*eta*chi**(25/9) + 3477817715529865344*eta**2*chi**(25/9) - 10967800378519109632*chi**(31/9) + 15859103224796779584*eta*chi**(31/9) - 570361367148627600*np.pi**2*eta*chi**(31/9) + 75826894562974464*eta**2*chi**(31/9) + 493039468942823701*chi**(38/9) - 6093524899450135056*eta*chi**(38/9) + 7224047870706816240*eta**2*chi**(38/9) + 8936345368795116909*chi**(44/9) - 13792229587322482296*eta*chi**(44/9) + 1646334720273668880*eta**2*chi**(44/9) + 9145059780731944160*chi**(50/9) - 15256596481285635840*eta*chi**(50/9) + 405971677516073760*np.pi**2*eta*chi**(50/9) - 109530449579765760*eta**2*chi**(50/9))) )
    
        #eccentricity
        x = []
        for n in [0,1,2,3,4,-4,-3,-2,-1]: 
            x.append( np.abs((G*M*2*np.pi*f)/( C**3*(l - (l + n)*k/(1 + k))))**(2/3) )
        x = np.array(x)

        et = (x**(3/2) * (et0**5 * ((94739615555*np.pi)/(958169088*chi**(113/18))-(1586634546601*np.pi)/(27786903552*chi**(95/18))-(1422200801*np.pi)/(13307904*chi**(25/6))+(1318556431*np.pi)/(26615808*chi**(19/6))+(607032981553*np.pi)/(27786903552*chi**(37/18))-(6029825087*np.pi)/(958169088*chi**(19/18)))+et0**3 * (-((1252771*np.pi)/(87552*chi**(25/6)))+(396797*np.pi)/(43776*chi**(19/6))+(1315151*np.pi)/(131328*chi**(37/18))-(1252771*np.pi)/(262656*chi**(19/18)))+et0 * ((377*np.pi)/(144*chi**(37/18))-(377*np.pi)/(144*chi**(19/18))))+\
    \
    x * (et0 * ((2833/2016-(197*eta)/72)/chi**(31/18)+(-(2833/2016)+(197*eta)/72)/chi**(19/18))+et0**3 * ((-(9414059/1225728)+(654631*eta)/43776)/chi**(23/6)+(11412055/5311488-(378697*eta)/43776)/chi**(19/6)+(386822573/47803392-(1482433*eta)/131328)/chi**(31/18)+(-(9414059/3677184)+(654631*eta)/131328)/chi**(19/18))+et0**5 * ((711929259595/13414367232-(49505846855*eta)/479084544)/chi**(107/18)+(-(1182747028465/174386774016)+(24493152461*eta)/479084544)/chi**(95/18)+(-(699589093187/9688154112)+(3092267495*eta)/26615808)/chi**(23/6)+(37922258765/3229384704-(1258410131*eta)/26615808)/chi**(19/6)+(3061519891285/174386774016-(11147601665*eta)/479084544)/chi**(31/18)+(-(45311656423/13414367232)+(3150863507*eta)/479084544)/chi**(19/18)))+\
        \
        x**(5/2) * (et0**5 * (((158961967498087*np.pi)/275952697344-(773693508027443*np.pi*eta)/482917220352)/chi**(125/18)+(-((171498319127425*np.pi)/1931668881408)+(46169592388985*np.pi*eta)/68988174336)/chi**(113/18)+(-((22474678352603165*np.pi)/56018397560832)+(1562835028401985*np.pi*eta)/2000657055744)/chi**(107/18)+((966616778184203239*np.pi)/2392785838669824-(102851838418889831*np.pi*eta)/322105785974784)/chi**(95/18)+(-((156129856813559*np.pi)/178858229760)+(70761454761383*np.pi*eta)/44714557440)/chi**(29/6)+((697744454755*np.pi)/5536088064-(162077392939*np.pi*eta)/319389696)/chi**(25/6)+((83537422031093*np.pi)/232515698688-(369245400305*np.pi*eta)/638779392)/chi**(23/6)+(-((6837459307169*np.pi)/19873136640)+(7877941103477*np.pi*eta)/44714557440)/chi**(19/6)+((3611887558130624419*np.pi)/11963929193349120-(517791739629486467*np.pi*eta)/1610528929873920)/chi**(49/18)+(-((1719724436739649*np.pi)/56018397560832)+(119585497365941*np.pi*eta)/2000657055744)/chi**(37/18)+(-((88784076847265*np.pi)/1931668881408)+(4202645827705*np.pi*eta)/68988174336)/chi**(31/18)+((158367949859977*np.pi)/9658344407040-(3240255264059*np.pi*eta)/2414586101760)/chi**(19/18))+et0**3 * ((-((12693032573*np.pi)/294174720)+(11292740311*np.pi*eta)/73543680)/chi**(29/6)+((330949595*np.pi)/19611648-(142768769*np.pi*eta)/2101248)/chi**(25/6)+((1124125901*np.pi)/29417472-(78169009*np.pi*eta)/1050624)/chi**(23/6)+(-((2057616403*np.pi)/32686080)+(2370731599*np.pi*eta)/73543680)/chi**(19/6)+((195499289159*np.pi)/2647572480-(65776041763*np.pi*eta)/661893120)/chi**(49/18)+(-((3725822783*np.pi)/264757248)+(259084747*np.pi*eta)/9455616)/chi**(37/18)+(-((11217854617*np.pi)/529514496)+(558877241*np.pi*eta)/18911232)/chi**(31/18)+((32902907141*np.pi)/2647572480-(673203247*np.pi*eta)/661893120)/chi**(19/18))+et0 * (((778843*np.pi)/1451520-(4996241*np.pi*eta)/362880)/chi**(49/18)+(-((1068041*np.pi)/290304)+(74269*np.pi*eta)/10368)/chi**(37/18)+(-((1068041*np.pi)/290304)+(74269*np.pi*eta)/10368)/chi**(31/18)+((9901567*np.pi)/1451520-(202589*np.pi*eta)/362880)/chi**(19/18)))+\
            \
            x**2 * (et0 * ((-(28850671/24385536)+(27565*eta)/145152+(33811*eta**2)/10368)/chi**(43/18)+(-(8025889/4064256)+(558101*eta)/72576-(38809*eta**2)/5184)/chi**(31/18)+(77006005/24385536-(1143767*eta)/145152+(43807*eta**2)/10368)/chi**(19/18))+et0**3 * ((-(9164199307/2118057984)+(1205846917*eta)/29417472-(13714021*eta**2)/233472)/chi**(9/2)+(32330351815/3569319936-(10345778159*eta)/191213568+(74603309*eta**2)/1050624)/chi**(23/6)+(8180980796033/1349202935808+(14604819923*eta)/2676989952-(317361763*eta**2)/14708736)/chi**(19/6)+(-(20952382669619/4047608807424)-(385200824731*eta)/24092909568+(4301644427*eta**2)/132378624)/chi**(43/18)+(-(1095868349309/96371638272)+(65400285919*eta)/1720922112-(292039301*eta**2)/9455616)/chi**(31/18)+(255890954615/44479217664-(3800737741*eta)/264757248+(145570661*eta**2)/18911232)/chi**(19/18))+et0**5 * ((16952610560003855/162260186038272-(79153315354555*eta)/137976348672+(47507268174605*eta**2)/68988174336)/chi**(119/18)+(-(16753611658206725/351563736416256)+(2837648691484435*eta)/6277923864576-(24125755174085*eta**2)/34494087168)/chi**(107/18)+(-(7937050519029473999/191953800083275776)-(1089957759112387*eta)/87890934104064+(1121044759543031*eta**2)/6277923864576)/chi**(95/18)+(-(25186092424407371/273438461657088)+(936816311138573*eta)/1627609890816-(1951606822255*eta**2)/2980970496)/chi**(9/2)+(2402572738143295/28211904774144-(55792908667709*eta)/116257849344+(352402173805*eta**2)/638779392)/chi**(23/6)+(27185399185217659/820315384971264+(48531816604129*eta)/1627609890816-(1054593138449*eta**2)/8942911488)/chi**(19/6)+(-(4719697288288984795/191953800083275776)-(4676818769915975*eta)/87890934104064+(669607180808035*eta**2)/6277923864576)/chi**(43/18)+(-(8673285852010405/351563736416256)+(506837220151715*eta)/6277923864576-(2196077528005*eta**2)/34494087168)/chi**(31/18)+(1231651832357155/162260186038272-(18293673608177*eta)/965834440704+(700659277417*eta**2)/68988174336)/chi**(19/18)))+\
                \
                et0**3 * (-(3323/(1824*chi**(19/6)))+3323/(1824*chi**(19/18)))+et0**5 * (50259743/(6653952*chi**(95/18))-11042329/(1108992*chi**(19/6))+15994231/(6653952*chi**(19/18)))+et0/chi**(19/18)+\
                \
                x**3 * (et0**3 * ((31472267987495/6167784849408-(318662569276073*eta)/4625838637056+(4844584781833*eta**2)/18356502528-(1562882519*eta**3)/5603328)/chi**(9/2)+(149592469*np.pi**2)/(2101248*chi**(25/6))+(23176718595161489/906664372862976-(866895029665039*eta)/32380870459392-(5814138473063*eta**2)/42831839232+(62520267311*eta**3)/353009664)/chi**(23/6)+(59358100103030627/8159979355766784+(2420024232862595*eta)/291427834134528-(103398129181999*eta**2)/1156459659264+(847423952119*eta**3)/9531260928)/chi**(43/18)-(495811927*np.pi**2)/(18911232*chi**(37/18))+(29787660990550865/1165711336538112-(591234360321013*eta)/5947506819072+(107636760191*eta**2)/874119168-(64940942431*eta**3)/1361608704)/chi**(31/18)+1/chi**(19/6)*(152896024020300184249/67999827964723200-(95207357*np.pi**2)/8404992-(245954159*gamma)/766080+(12374839994637661*eta)/10793623486464-(116237911*np.pi**2*eta)/1400832-(3908281091711*eta**2)/128495517696-(42680326813*eta**3)/1059028992-(33962745773*np.log(2))/2298240+(5362264233*np.log(3))/680960-(245954159*np.log(x))/1532160)+1/chi**(19/18)*(-(110724557880778937/704550807797760)+(600535883*np.pi**2)/75644928+(11022391*gamma)/459648+(536131194179051*eta)/16012518359040+(13215571*np.pi**2*eta)/4202496-(1193082406697*eta**2)/38125043712+(35382609493*eta**3)/4084826112+(40178393*np.log(2))/6894720+(28800441*np.log(3))/680960+(11022391*np.log(x))/919296)+1/chi**(55/18)*(-(3881667007528080426037/2243994322835865600)+(720177509*np.pi**2)/75644928+(517414657*gamma)/2298240-(1395931720786001359*eta)/1457139170672640+(295851449*np.pi**2*eta)/4202496-(112681906698415*eta**2)/3469378977792-(1549239851389*eta**3)/28593782784+(101727523747*np.log(2))/6894720-(5477465997*np.log(3))/680960+(517414657*np.log(x))/4596480-(517414657*np.log(chi))/6894720)+1/chi**(31/6)*(-(99813874374700537/234850269265920)-(429547595*np.pi**2)/8404992+(11022391*gamma)/153216-(62659748948903*eta)/1779168706560+(13215571*np.pi**2*eta)/1400832-(95613034561*eta**2)/1412038656+(22151672941*eta**3)/151289856+(40178393*np.log(2))/2298240+(86401323*np.log(3))/680960+(11022391*np.log(x))/306432-(11022391*np.log(chi))/459648))+et0 * ((81733950943/49161240576-(6152132057*eta)/1755758592-(1348031*eta**2)/331776+(6660767*eta**3)/746496)/chi**(43/18)-(142129*np.pi**2)/(20736*chi**(37/18))+(218158012165/49161240576-(34611934451*eta)/1755758592+(191583143*eta**2)/6967296-(8629979*eta**3)/746496)/chi**(31/18)+1/chi**(19/18)*(-(33320661414619/386266890240)+(180721*np.pi**2)/41472+(3317*gamma)/252+(161339510737*eta)/8778792960+(3977*np.pi**2*eta)/2304-(359037739*eta**2)/20901888+(10647791*eta**3)/2239488+(12091*np.log(2))/3780+(26001*np.log(3))/1120+(3317*np.log(x))/504)+1/chi**(55/18)*(216750571931393/2703868231680+(103537*np.pi**2)/41472-(3317*gamma)/252+(866955547*eta)/179159040-(3977*np.pi**2*eta)/2304-(130785737*eta**2)/20901888-(4740155*eta**3)/2239488-(12091*np.log(2))/3780-(26001*np.log(3))/1120-(3317*np.log(x))/504+(3317*np.log(chi))/756))+et0**5 * ((-(398940554960039073025/4252514955691032576)+(751550699250552365*eta)/614880704987136-(2712807799571918515*eta**2)/602680690999296+(23151784966473335*eta**3)/4967148552192)/chi**(119/18)-(103131245529065*np.pi**2)/(137976348672*chi**(113/18))+(-(112428320602052499195835/386978860967883964416)+(6613733131933528864325*eta)/13820673605995855872+(6008938601459478835*eta**2)/4218764836995072-(1104229088149885535*eta**3)/452010518249472)/chi**(107/18)+(86495658134944405735/796252800345440256-(5011810129939409773*eta)/4490147370369024+(273844614114965969*eta**2)/78125274759168-(222409764901445*eta**3)/71543291904)/chi**(9/2)+(169823957639*np.pi**2)/(319389696*chi**(25/6))+(1722336724789945193177/7166275203108962304-(43622763718525097819*eta)/255938400111034368-(93995440212410537*eta**2)/78125274759168+(295325748986095*eta**3)/214629875712)/chi**(23/6)+(13370902417722693924235/386978860967883964416+(5456578161604350265*eta)/727403873999781888-(95888813809642495*eta**2)/324520372076544+(131912614619182895*eta**3)/452010518249472)/chi**(43/18)-(7891428760189*np.pi**2)/(137976348672*chi**(37/18))+(235755416055892166425/4252514955691032576-(4595658861880171685*eta)/21696504875974656+(22172511125256925*eta**2)/86097241571328-(488342986138655*eta**3)/4967148552192)/chi**(31/18)+1/chi**(19/6)*(508073487819457512259427/41343895402551705600-(316374047311*np.pi**2)/5110235136-(817305670357*gamma)/465776640+(41121593302180947503*eta)/6562523079770112-(386258578253*np.pi**2*eta)/851705856-(12987218067755653*eta**2)/78125274759168-(141826725999599*eta**3)/643889627136-(112858204203679*np.log(2))/1397329920+(17818804046259*np.log(3))/414023680-(817305670357*np.log(x))/931553280)+1/chi**(19/18)*(-(48448941430745732999/233654667895111680)+(2890493420551*np.pi**2)/275952697344+(53052864227*gamma)/1676795904+(2580501404154558247*eta)/58413666973777920+(63609056687*np.pi**2*eta)/15330705408-(5742532535283709*eta**2)/139080159461376+(170303228893721*eta**3)/14901445656576+(193386247021*np.log(2))/25151938560+(46207333359*np.log(3))/828047360+(53052864227*np.log(x))/3353591808)+1/chi**(95/18)*(-(1381494982108367906597498839/106419186766168090214400)+(210263419125757*np.pi**2)/1379763486720+(76993944487871*gamma)/41919897600-(459744548596452548880253*eta)/69103368029979279360+(36313178438959*np.pi**2*eta)/76653527040+(15325478660942336659*eta**2)/63281472554926080+(602052690061119545*eta**3)/1356031554748416+(11797275511817693*np.log(2))/125759692800-(196321353899121*np.log(3))/4140236800-(208984375*np.log(5))/86016+(76993944487871*np.log(x))/83839795200)+1/chi**(55/18)*(-(1077842077038057332046031457/106419186766168090214400)+(22839987538133*np.pi**2)/1379763486720+(56669035748389*gamma)/41919897600-(332702882500991633514587*eta)/69103368029979279360+(29987633672981*np.pi**2*eta)/76653527040-(9909009804994157899*eta**2)/63281472554926080-(394569227909087057*eta**3)/1356031554748416+(8507394273435727*np.log(2))/125759692800-(162397985681259*np.log(3))/4140236800+(208984375*np.log(5))/86016+(56669035748389*np.log(x))/83839795200-(56669035748389*np.log(chi))/125759692800)+1/chi**(131/18)*(10744314227409210983279/3598281885584719872+(168886055311895*np.pi**2)/275952697344-(833557837655*gamma)/1676795904-(854950976221006253*eta)/1668961913536512-(999414989555*np.pi**2*eta)/15330705408+(384051995415241885*eta**2)/139080159461376-(44332570043722025*eta**3)/14901445656576-(607690552613*np.log(2))/5030387712-(145200397527*np.log(3))/165609472-(833557837655*np.log(x))/3353591808+(833557837655*np.log(chi))/5030387712)+1/chi**(31/6)*(3632678520819467743258727/454782849428068761600-(2333251524667*np.pi**2)/5110235136-(451031617427*gamma)/465776640+(41341286137777666139*eta)/7572142015119360-(298427893387*np.pi**2*eta)/851705856-(142707649066100323*eta**2)/78125274759168+(1514475909193295*eta**3)/643889627136-(112591178603801*np.log(2))/1397329920+(18393027238917*np.log(3))/414023680-(451031617427*np.log(x))/931553280+(451031617427*np.log(chi))/1397329920))))
        
        #NOTE: et(n value, l-1)
        #I could have used a transpose but I am saving computation time

        ###################cplus0_start####################
        #Cp0[l-1,n]
        Cp0=np.zeros((10,3))
    
        Cp0[0,-2]=et[0,-2]**3*(-(13/16)*cosine(2*beta)-13/16*cosine(iota)**2*cosine(2*beta))+et[0,-2]**5*(-(5/384)*cosine(2*beta)-5/384*cosine(iota)**2*cosine(2*beta))+et[0,-2]*(3/2*cosine(2*beta)+3/2*cosine(iota)**2*cosine(2*beta))

        Cp0[1,-2]=-2*cosine(2*beta)-2*cosine(iota)**2*cosine(2*beta)+et[1,-2]**4*(-(23/8)*cosine(2*beta)-23/8*cosine(iota)**2*cosine(2*beta))+et[1,-2]**6*(65/144*cosine(2*beta)+65/144*cosine(iota)**2*cosine(2*beta))+et[1,-2]**2*(5*cosine(2*beta)+5*cosine(iota)**2*cosine(2*beta))

        Cp0[2,-2]=et[2,-2]**5*(-(963/128)*cosine(2*beta)-963/128*cosine(iota)**2*cosine(2*beta))+et[2,-2]*(-(9/2)*cosine(2*beta)-9/2*cosine(iota)**2*cosine(2*beta))+et[2,-2]**3*(171/16*cosine(2*beta)+171/16*cosine(iota)**2*cosine(2*beta))

        Cp0[3,-2]=et[3,-2]**6*(-(101/6)*cosine(2*beta)-101/6*cosine(iota)**2*cosine(2*beta))+et[3,-2]**2*(-8*cosine(2*beta)-8*cosine(iota)**2*cosine(2*beta))+et[3,-2]**4*(20*cosine(2*beta)+20*cosine(iota)**2*cosine(2*beta))

        Cp0[4,-2]=et[4,-2]**3*(-(625/48)*cosine(2*beta)-625/48*cosine(iota)**2*cosine(2*beta))+et[4,-2]**5*(26875/768*cosine(2*beta)+26875/768*cosine(iota)**2*cosine(2*beta))

        Cp0[5,-2]=et[5,-2]**4*(-(81/4)*cosine(2*beta)-81/4*cosine(iota)**2*cosine(2*beta))+et[5,-2]**6*(2349/40*cosine(2*beta)+2349/40*cosine(iota)**2*cosine(2*beta))

        Cp0[6,-2]=et[6,-2]**5*(-((117649*cosine(2*beta))/3840)-(117649*cosine(iota)**2*cosine(2*beta))/3840)

        Cp0[7,-2]=et[7,-2]**6*(-(2048/45)*cosine(2*beta)-2048/45*cosine(iota)**2*cosine(2*beta))


        Cp0[0,0]=et[0,0]*sine(iota)**2-1/8*et[0,0]**3*sine(iota)**2+1/192*et[0,0]**5*sine(iota)**2

        Cp0[1,0]=et[1,0]**2*sine(iota)**2-1/3*et[1,0]**4*sine(iota)**2+1/24*et[1,0]**6*sine(iota)**2

        Cp0[2,0]=9/8*et[2,0]**3*sine(iota)**2-81/128*et[2,0]**5*sine(iota)**2

        Cp0[3,0]=4/3*et[3,0]**4*sine(iota)**2-16/15*et[3,0]**6*sine(iota)**2

        Cp0[4,0]=625/384*et[4,0]**5*sine(iota)**2

        Cp0[5,0]=81/40*et[5,0]**6*sine(iota)**2


        Cp0[0,2]=et[0,2]**5*(47/768*cosine(2*beta)+47/768*cosine(iota)**2*cosine(2*beta))+et[0,2]**3*(7/48*cosine(2*beta)+7/48*cosine(iota)**2*cosine(2*beta))

        Cp0[1,2]=et[1,2]**6*(11/240*cosine(2*beta)+11/240*cosine(iota)**2*cosine(2*beta))+et[1,2]**4*(1/8*cosine(2*beta)+1/8*cosine(iota)**2*cosine(2*beta))

        Cp0[2,2]=et[2,2]**5*((153*cosine(2*beta))/1280+(153*cosine(iota)**2*cosine(2*beta))/1280)

        Cp0[3,2]=et[3,2]**6*(11/90*cosine(2*beta)+11/90*cosine(iota)**2*cosine(2*beta))
        
        ###################splus0_start####################
        Sp0=np.zeros((10,3))

        Sp0[0,-2]=et[0,-2]**3*(-(13/16)*sine(2*beta)-13/16*cosine(iota)**2*sine(2*beta))+et[0,-2]**5*(-(5/384)*sine(2*beta)-5/384*cosine(iota)**2*sine(2*beta))+et[0,-2]*(3/2*sine(2*beta)+3/2*cosine(iota)**2*sine(2*beta))

        Sp0[1,-2]=-2*sine(2*beta)-2*cosine(iota)**2*sine(2*beta)+et[1,-2]**4*(-(23/8)*sine(2*beta)-23/8*cosine(iota)**2*sine(2*beta))+et[1,-2]**6*(65/144*sine(2*beta)+65/144*cosine(iota)**2*sine(2*beta))+et[1,-2]**2*(5*sine(2*beta)+5*cosine(iota)**2*sine(2*beta))

        Sp0[2,-2]=et[2,-2]**5*(-(963/128)*sine(2*beta)-963/128*cosine(iota)**2*sine(2*beta))+et[2,-2]*(-(9/2)*sine(2*beta)-9/2*cosine(iota)**2*sine(2*beta))+et[2,-2]**3*(171/16*sine(2*beta)+171/16*cosine(iota)**2*sine(2*beta))

        Sp0[3,-2]=et[3,-2]**6*(-(101/6)*sine(2*beta)-101/6*cosine(iota)**2*sine(2*beta))+et[3,-2]**2*(-8*sine(2*beta)-8*cosine(iota)**2*sine(2*beta))+et[3,-2]**4*(20*sine(2*beta)+20*cosine(iota)**2*sine(2*beta))

        Sp0[4,-2]=et[4,-2]**3*(-(625/48)*sine(2*beta)-625/48*cosine(iota)**2*sine(2*beta))+et[4,-2]**5*(26875/768*sine(2*beta)+26875/768*cosine(iota)**2*sine(2*beta))

        Sp0[5,-2]=et[5,-2]**4*(-(81/4)*sine(2*beta)-81/4*cosine(iota)**2*sine(2*beta))+et[5,-2]**6*(2349/40*sine(2*beta)+2349/40*cosine(iota)**2*sine(2*beta))

        Sp0[6,-2]=et[6,-2]**5*(-((117649*sine(2*beta))/3840)-(117649*cosine(iota)**2*sine(2*beta))/3840)

        Sp0[7,-2]=et[7,-2]**6*(-(2048/45)*sine(2*beta)-2048/45*cosine(iota)**2*sine(2*beta))


        Sp0[0,2]=et[0,2]**3*(-(7/48)*sine(2*beta)-7/48*cosine(iota)**2*sine(2*beta))+et[0,2]**5*(-(47/768)*sine(2*beta)-47/768*cosine(iota)**2*sine(2*beta))

        Sp0[1,2]=et[1,2]**4*(-(1/8)*sine(2*beta)-1/8*cosine(iota)**2*sine(2*beta))+et[1,2]**6*(-(11/240)*sine(2*beta)-11/240*cosine(iota)**2*sine(2*beta))

        Sp0[2,2]=et[2,2]**5*(-((153*sine(2*beta))/1280)-(153*cosine(iota)**2*sine(2*beta))/1280)

        Sp0[3,2]=et[3,2]**6*(-(11/90)*sine(2*beta)-11/90*cosine(iota)**2*sine(2*beta))
        
        ###################ccross0_start####################
        Cx0=np.zeros((10,3))

        Cx0[0,-2]=-3*et[0,-2]*cosine(iota)*sine(2*beta)+13/8*et[0,-2]**3*cosine(iota)*sine(2*beta)+5/192*et[0,-2]**5*cosine(iota)*sine(2*beta)

        Cx0[1,-2]=4*cosine(iota)*sine(2*beta)-10*et[1,-2]**2*cosine(iota)*sine(2*beta)+23/4*et[1,-2]**4*cosine(iota)*sine(2*beta)-65/72*et[1,-2]**6*cosine(iota)*sine(2*beta)

        Cx0[2,-2]=9*et[2,-2]*cosine(iota)*sine(2*beta)-171/8*et[2,-2]**3*cosine(iota)*sine(2*beta)+963/64*et[2,-2]**5*cosine(iota)*sine(2*beta)

        Cx0[3,-2]=16*et[3,-2]**2*cosine(iota)*sine(2*beta)-40*et[3,-2]**4*cosine(iota)*sine(2*beta)+101/3*et[3,-2]**6*cosine(iota)*sine(2*beta)

        Cx0[4,-2]=625/24*et[4,-2]**3*cosine(iota)*sine(2*beta)-26875/384*et[4,-2]**5*cosine(iota)*sine(2*beta)

        Cx0[5,-2]=81/2*et[5,-2]**4*cosine(iota)*sine(2*beta)-2349/20*et[5,-2]**6*cosine(iota)*sine(2*beta)

        Cx0[6,-2]=(117649*et[6,-2]**5*cosine(iota)*sine(2*beta))/1920

        Cx0[7,-2]=4096/45*et[7,-2]**6*cosine(iota)*sine(2*beta)


        Cx0[0,2]=-(7/24)*et[0,2]**3*cosine(iota)*sine(2*beta)-47/384*et[0,2]**5*cosine(iota)*sine(2*beta)

        Cx0[1,2]=-(1/4)*et[1,2]**4*cosine(iota)*sine(2*beta)-11/120*et[1,2]**6*cosine(iota)*sine(2*beta)

        Cx0[2,2]=-(153/640)*et[2,2]**5*cosine(iota)*sine(2*beta)

        Cx0[3,2]=-(11/45)*et[3,2]**6*cosine(iota)*sine(2*beta)
        
        ###################scross0_start####################
        Sx0=np.zeros((10,3))
        
        Sx0[0,-2]=3*et[0,-2]*cosine(iota)*cosine(2*beta)-13/8*et[0,-2]**3*cosine(iota)*cosine(2*beta)-5/192*et[0,-2]**5*cosine(iota)*cosine(2*beta)

        Sx0[1,-2]=-4*cosine(iota)*cosine(2*beta)+10*et[1,-2]**2*cosine(iota)*cosine(2*beta)-23/4*et[1,-2]**4*cosine(iota)*cosine(2*beta)+65/72*et[1,-2]**6*cosine(iota)*cosine(2*beta)

        Sx0[2,-2]=-9*et[2,-2]*cosine(iota)*cosine(2*beta)+171/8*et[2,-2]**3*cosine(iota)*cosine(2*beta)-963/64*et[2,-2]**5*cosine(iota)*cosine(2*beta)

        Sx0[3,-2]=-16*et[3,-2]**2*cosine(iota)*cosine(2*beta)+40*et[3,-2]**4*cosine(iota)*cosine(2*beta)-101/3*et[3,-2]**6*cosine(iota)*cosine(2*beta)

        Sx0[4,-2]=-(625/24)*et[4,-2]**3*cosine(iota)*cosine(2*beta)+26875/384*et[4,-2]**5*cosine(iota)*cosine(2*beta)

        Sx0[5,-2]=-(81/2)*et[5,-2]**4*cosine(iota)*cosine(2*beta)+2349/20*et[5,-2]**6*cosine(iota)*cosine(2*beta)

        Sx0[6,-2]=-((117649*et[6,-2]**5*cosine(iota)*cosine(2*beta))/1920)

        Sx0[7,-2]=-(4096/45)*et[7,-2]**6*cosine(iota)*cosine(2*beta)


        Sx0[0,2]=-(7/24)*et[0,2]**3*cosine(iota)*cosine(2*beta)-47/384*et[0,2]**5*cosine(iota)*cosine(2*beta)

        Sx0[1,2]=-(1/4)*et[1,2]**4*cosine(iota)*cosine(2*beta)-11/120*et[1,2]**6*cosine(iota)*cosine(2*beta)

        Sx0[2,2]=-(153/640)*et[2,2]**5*cosine(iota)*cosine(2*beta)

        Sx0[3,2]=-(11/45)*et[3,2]**6*cosine(iota)*cosine(2*beta)
        

        ###################cplus05_start####################
        Cp05 = np.zeros((10,5))

        Cp05[0,1]=et[0,1]**6*((109*cosine(beta)*sine(iota))/12288+(265*cosine(iota)**2*cosine(beta)*sine(iota))/12288)+et[0,1]**4*(1/64*cosine(beta)*sine(iota)+7/192*cosine(iota)**2*cosine(beta)*sine(iota))+et[0,1]**2*(-(9/32)*cosine(beta)*sine(iota)+11/32*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[1,1]=et[1,1]**5*(1/16*cosine(beta)*sine(iota)-1/48*cosine(iota)**2*cosine(beta)*sine(iota))+et[1,1]**3*(-(1/3)*cosine(beta)*sine(iota)+1/3*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[2,1]=et[2,1]**6*((369*cosine(beta)*sine(iota))/2560-(243*cosine(iota)**2*cosine(beta)*sine(iota))/2560)+et[2,1]**4*(-(207/512)*cosine(beta)*sine(iota)+189/512*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[3,1]=et[3,1]**5*(-(1/2)*cosine(beta)*sine(iota)+13/30*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[4,1]=et[4,1]**6*(-((23125*cosine(beta)*sine(iota))/36864)+(19375*cosine(iota)**2*cosine(beta)*sine(iota))/36864)


        Cp05[0,-1]=-(5/4)*cosine(beta)*sine(iota)-1/4*cosine(iota)**2*cosine(beta)*sine(iota)+et[0,-1]**2*(3/2*cosine(beta)*sine(iota)-1/2*cosine(iota)**2*cosine(beta)*sine(iota))+et[0,-1]**6*((41*cosine(beta)*sine(iota))/2304+(37*cosine(iota)**2*cosine(beta)*sine(iota))/2304)+et[0,-1]**4*(-(51/256)*cosine(beta)*sine(iota)+41/256*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[1,-1]=et[1,-1]**3*(4*cosine(beta)*sine(iota)-2*cosine(iota)**2*cosine(beta)*sine(iota))+et[1,-1]**5*(-(19/16)*cosine(beta)*sine(iota)+35/48*cosine(iota)**2*cosine(beta)*sine(iota))+et[1,-1]*(-3*cosine(beta)*sine(iota)+cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[2,-1]=et[2,-1]**4*(531/64*cosine(beta)*sine(iota)-297/64*cosine(iota)**2*cosine(beta)*sine(iota))+et[2,-1]**6*(-((15399*cosine(beta)*sine(iota))/4096)+(9477*cosine(iota)**2*cosine(beta)*sine(iota))/4096)+et[2,-1]**2*(-(171/32)*cosine(beta)*sine(iota)+81/32*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[3,-1]=et[3,-1]**5*(31/2*cosine(beta)*sine(iota)-55/6*cosine(iota)**2*cosine(beta)*sine(iota))+et[3,-1]**3*(-(26/3)*cosine(beta)*sine(iota)+14/3*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[4,-1]=et[4,-1]**6*((41875*cosine(beta)*sine(iota))/1536-(25625*cosine(iota)**2*cosine(beta)*sine(iota))/1536)+et[4,-1]**4*(-(6875/512)*cosine(beta)*sine(iota)+(11875*cosine(iota)**2*cosine(beta)*sine(iota))/1536)
        Cp05[5,-1]=et[5,-1]**5*(-(81/4)*cosine(beta)*sine(iota)+243/20*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp05[6,-1]=et[6,-1]**6*(-((5529503*cosine(beta)*sine(iota))/184320)+(3411821*cosine(iota)**2*cosine(beta)*sine(iota))/184320)
        
        
        Cp05[0,3]=et[0,3]**4*(-(25/512)*cosine(3*beta)*sine(iota)-25/512*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[0,3]**6*(-((179*cosine(3*beta)*sine(iota))/7680)-(179*cosine(iota)**2*cosine(3*beta)*sine(iota))/7680)
        Cp05[1,3]=et[1,3]**5*(-(13/240)*cosine(3*beta)*sine(iota)-13/240*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[2,3]=et[2,3]**6*(-((1233*cosine(3*beta)*sine(iota))/20480)-(1233*cosine(iota)**2*cosine(3*beta)*sine(iota))/20480)
      
    
        Cp05[0,-3]=et[0,-3]**4*(-(65/192)*cosine(3*beta)*sine(iota)-65/192*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[0,-3]**6*(-((19*cosine(3*beta)*sine(iota))/12288)-(19*cosine(iota)**2*cosine(3*beta)*sine(iota))/12288)+et[0,-3]**2*(19/32*cosine(3*beta)*sine(iota)+19/32*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[1,-3]=et[1,-3]*(-3*cosine(3*beta)*sine(iota)-3*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[1,-3]**5*(-(133/48)*cosine(3*beta)*sine(iota)-133/48*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[1,-3]**3*(11/2*cosine(3*beta)*sine(iota)+11/2*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[2,-3]=9/4*cosine(3*beta)*sine(iota)+9/4*cosine(iota)**2*cosine(3*beta)*sine(iota)+et[2,-3]**2*(-(27/2)*cosine(3*beta)*sine(iota)-27/2*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[2,-3]**6*(-(3141/256)*cosine(3*beta)*sine(iota)-3141/256*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[2,-3]**4*(5319/256*cosine(3*beta)*sine(iota)+5319/256*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[3,-3]=et[3,-3]**3*(-38*cosine(3*beta)*sine(iota)-38*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[3,-3]*(8*cosine(3*beta)*sine(iota)+8*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[3,-3]**5*(176/3*cosine(3*beta)*sine(iota)+176/3*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[4,-3]=et[4,-3]**4*(-(5625/64)*cosine(3*beta)*sine(iota)-5625/64*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[4,-3]**2*(625/32*cosine(3*beta)*sine(iota)+625/32*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[4,-3]**6*((1746875*cosine(3*beta)*sine(iota))/12288+(1746875*cosine(iota)**2*cosine(3*beta)*sine(iota))/12288)
        Cp05[5,-3]=et[5,-3]**5*(-(729/4)*cosine(3*beta)*sine(iota)-729/4*cosine(iota)**2*cosine(3*beta)*sine(iota))+et[5,-3]**3*(81/2*cosine(3*beta)*sine(iota)+81/2*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[6,-3]=et[6,-3]**6*(-((2705927*cosine(3*beta)*sine(iota))/7680)-(2705927*cosine(iota)**2*cosine(3*beta)*sine(iota))/7680)+et[6,-3]**4*((117649*cosine(3*beta)*sine(iota))/1536+(117649*cosine(iota)**2*cosine(3*beta)*sine(iota))/1536)
        Cp05[7,-3]=et[7,-3]**5*(2048/15*cosine(3*beta)*sine(iota)+2048/15*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp05[8,-3]=et[8,-3]**6*((4782969*cosine(3*beta)*sine(iota))/20480+(4782969*cosine(iota)**2*cosine(3*beta)*sine(iota))/20480)  

        ###################splus05_start####################
        Sp05 = np.zeros((10,5))

        Sp05[0,-1] =  - (5/4)*sine(iota)*sine(beta) - 1/4*cosine(iota)**2*sine(iota)*sine(beta) + et[0,-1]**2*(3/2*sine(iota)*sine(beta) - 1/2*cosine(iota)**2*sine(iota)*sine(beta)) + et[0,-1]**6*((41*sine(iota)*sine(beta))/2304 + (37*cosine(iota)**2*sine(iota)*sine(beta))/2304) + et[0,-1]**4*( - (51/256)*sine(iota)*sine(beta) + 41/256*cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[1,-1] = et[1,-1]**3*(4*sine(iota)*sine(beta) - 2*cosine(iota)**2*sine(iota)*sine(beta)) + et[1,-1]**5*( - (19/16)*sine(iota)*sine(beta) + 35/48*cosine(iota)**2*sine(iota)*sine(beta)) + et[1,-1]*( - 3*sine(iota)*sine(beta) + cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[2,-1] = et[2,-1]**4*(531/64*sine(iota)*sine(beta) - 297/64*cosine(iota)**2*sine(iota)*sine(beta)) + et[2,-1]**6*( - ((15399*sine(iota)*sine(beta))/4096) + (9477*cosine(iota)**2*sine(iota)*sine(beta))/4096) + et[2,-1]**2*( - (171/32)*sine(iota)*sine(beta) + 81/32*cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[3,-1] = et[3,-1]**5*(31/2*sine(iota)*sine(beta) - 55/6*cosine(iota)**2*sine(iota)*sine(beta)) + et[3,-1]**3*( - (26/3)*sine(iota)*sine(beta) + 14/3*cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[4,-1] = et[4,-1]**6*((41875*sine(iota)*sine(beta))/1536 - (25625*cosine(iota)**2*sine(iota)*sine(beta))/1536) + et[4,-1]**4*( - (6875/512)*sine(iota)*sine(beta) + (11875*cosine(iota)**2*sine(iota)*sine(beta))/1536)

        Sp05[5,-1] = et[5,-1]**5*( - (81/4)*sine(iota)*sine(beta) + 243/20*cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[6,-1] = et[6,-1]**6*( - ((5529503*sine(iota)*sine(beta))/184320) + (3411821*cosine(iota)**2*sine(iota)*sine(beta))/184320)


        Sp05[0,1] = et[0,1]**2*(9/32*sine(iota)*sine(beta) - 11/32*cosine(iota)**2*sine(iota)*sine(beta)) + et[0,1]**4*( - (1/64)*sine(iota)*sine(beta) - 7/192*cosine(iota)**2*sine(iota)*sine(beta)) + et[0,1]**6*( - ((109*sine(iota)*sine(beta))/12288) - (265*cosine(iota)**2*sine(iota)*sine(beta))/12288)

        Sp05[1,1] = et[1,1]**3*(1/3*sine(iota)*sine(beta) - 1/3*cosine(iota)**2*sine(iota)*sine(beta)) + et[1,1]**5*( - (1/16)*sine(iota)*sine(beta) + 1/48*cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[2,1] = et[2,1]**4*(207/512*sine(iota)*sine(beta) - 189/512*cosine(iota)**2*sine(iota)*sine(beta)) + et[2,1]**6*( - ((369*sine(iota)*sine(beta))/2560) + (243*cosine(iota)**2*sine(iota)*sine(beta))/2560)

        Sp05[3,1] = et[3,1]**5*(1/2*sine(iota)*sine(beta) - 13/30*cosine(iota)**2*sine(iota)*sine(beta))

        Sp05[4,1] = et[4,1]**6*((23125*sine(iota)*sine(beta))/36864 - (19375*cosine(iota)**2*sine(iota)*sine(beta))/36864)


        Sp05[0,-3] = et[0,-3]**4*( - (65/192)*sine(iota)*sine(3*beta) - 65/192*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[0,-3]**6*( - ((19*sine(iota)*sine(3*beta))/12288) - (19*cosine(iota)**2*sine(iota)*sine(3*beta))/12288) + et[0,-3]**2*(19/32*sine(iota)*sine(3*beta) + 19/32*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[1,-3] = et[1,-3]*( - 3*sine(iota)*sine(3*beta) - 3*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[1,-3]**5*( - (133/48)*sine(iota)*sine(3*beta) - 133/48*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[1,-3]**3*(11/2*sine(iota)*sine(3*beta) + 11/2*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[2,-3] = 9/4*sine(iota)*sine(3*beta) + 9/4*cosine(iota)**2*sine(iota)*sine(3*beta) + et[2,-3]**2*( - (27/2)*sine(iota)*sine(3*beta) - 27/2*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[2,-3]**6*( - (3141/256)*sine(iota)*sine(3*beta) - 3141/256*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[2,-3]**4*(5319/256*sine(iota)*sine(3*beta) + 5319/256*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[3,-3] = et[3,-3]**3*( - 38*sine(iota)*sine(3*beta) - 38*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[3,-3]**(8*sine(iota)*sine(3*beta) + 8*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[3,-3]**5*(176/3*sine(iota)*sine(3*beta) + 176/3*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[4,-3] = et[4,-3]**4*( - (5625/64)*sine(iota)*sine(3*beta) - 5625/64*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[4,-3]**2*(625/32*sine(iota)*sine(3*beta) + 625/32*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[4,-3]**6*((1746875*sine(iota)*sine(3*beta))/12288 + (1746875*cosine(iota)**2*sine(iota)*sine(3*beta))/12288)

        Sp05[5,-3] = et[5,-3]**5*( - (729/4)*sine(iota)*sine(3*beta) - 729/4*cosine(iota)**2*sine(iota)*sine(3*beta)) + et[5,-3]**3*(81/2*sine(iota)*sine(3*beta) + 81/2*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[6,-3] = et[6,-3]**6*( - ((2705927*sine(iota)*sine(3*beta))/7680) - (2705927*cosine(iota)**2*sine(iota)*sine(3*beta))/7680) + et[6,-3]**4*((117649*sine(iota)*sine(3*beta))/1536 + (117649*cosine(iota)**2*sine(iota)*sine(3*beta))/1536)

        Sp05[7,-3] = et[7,-3]**5*(2048/15*sine(iota)*sine(3*beta) + 2048/15*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[8,-3] = et[8,-3]**6*((4782969*sine(iota)*sine(3*beta))/20480 + (4782969*cosine(iota)**2*sine(iota)*sine(3*beta))/20480)


        Sp05[0,3] = et[0,3]**6*((179*sine(iota)*sine(3*beta))/7680 + (179*cosine(iota)**2*sine(iota)*sine(3*beta))/7680) + et[0,3]**4*(25/512*sine(iota)*sine(3*beta) + 25/512*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[1,3] = et[1,3]**5*(13/240*sine(iota)*sine(3*beta) + 13/240*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp05[2,3] = et[2,3]**6*((1233*sine(iota)*sine(3*beta))/20480 + (1233*cosine(iota)**2*sine(iota)*sine(3*beta))/20480)

        ###################ccross05_start####################
        Cx05 = np.zeros((10,5))

        Cx05[0,1]=-(1/16)*et[0,1]**2 *cosine(iota)*sine(iota)*sine(beta)-5/96*et[0,1]**4 *cosine(iota)*sine(iota)*sine(beta)-(187*et[0,1]**6 *cosine(iota)*sine(iota)*sine(beta))/6144
        Cx05[1,1]=-(1/24)*et[1,1]**5 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[2,1]=9/256*et[2,1]**4 *cosine(iota)*sine(iota)*sine(beta)-(63*et[2,1]**6 *cosine(iota)*sine(iota)*sine(beta))/1280
        Cx05[3,1]=1/15*et[3,1]**5 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[4,1]=(625*et[4,1]**6 *cosine(iota)*sine(iota)*sine(beta))/6144
        
        
        Cx05[0,-1]=3/2*cosine(iota)*sine(iota)*sine(beta)-et[0,-1]**2 *cosine(iota)*sine(iota)*sine(beta)+5/128*et[0,-1]**4 *cosine(iota)*sine(iota)*sine(beta)-13/384*et[0,-1]**6 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[1,-1]=2*et[1,-1]*cosine(iota)*sine(iota)*sine(beta)-2*et[1,-1]**3 *cosine(iota)*sine(iota)*sine(beta)+11/24*et[1,-1]**5 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[2,-1]=45/16*et[2,-1]**2 *cosine(iota)*sine(iota)*sine(beta)-117/32*et[2,-1]**4 *cosine(iota)*sine(iota)*sine(beta)+(2961*et[2,-1]**6 *cosine(iota)*sine(iota)*sine(beta))/2048
        Cx05[3,-1]=4*et[3,-1]**3 *cosine(iota)*sine(iota)*sine(beta)-19/3*et[3,-1]**5 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[4,-1]=4375/768*et[4,-1]**4 *cosine(iota)*sine(iota)*sine(beta)-8125/768*et[4,-1]**6 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[5,-1]=81/10*et[5,-1]**5 *cosine(iota)*sine(iota)*sine(beta)
        Cx05[6,-1]=(117649*et[6,-1]**6 *cosine(iota)*sine(iota)*sine(beta))/10240 
        
        
        Cx05[0,3]=25/256*et[0,3]**4 *cosine(iota)*sine(iota)*sine(3*beta)+(179*et[0,3]**6 *cosine(iota)*sine(iota)*sine(3*beta))/3840
        Cx05[1,3]=13/120*et[1,3]**5 *cosine(iota)*sine(iota)*sine(3*beta)
        Cx05[2,3]=(1233*et[2,3]**6 *cosine(iota)*sine(iota)*sine(3*beta))/10240        
        
        
        Cx05[0,-3]=-(19/16)*et[0,-3]**2 *cosine(iota)*sine(iota)*sine(3*beta)+65/96*et[0,-3]**4 *cosine(iota)*sine(iota)*sine(3*beta)+(19*et[0,-3]**6 *cosine(iota)*sine(iota)*sine(3*beta))/6144
        Cx05[1,-3]=6*et[1,-3]*cosine(iota)*sine(iota)*sine(3*beta)-11*et[1,-3]**3 *cosine(iota)*sine(iota)*sine(3*beta)+133/24*et[1,-3]**5 *cosine(iota)*sine(iota)*sine(3*beta)
        Cx05[2,-3]=-(9/2)*cosine(iota)*sine(iota)*sine(3*beta)+27*et[2,-3]**2 *cosine(iota)*sine(iota)*sine(3*beta)-5319/128*et[2,-3]**4 *cosine(iota)*sine(iota)*sine(3*beta)+3141/128*et[2,-3]**6 *cosine(iota)*sine(iota)*sine(3*beta)
        Cx05[3,-3]=-16*et[3,-3]*cosine(iota)*sine(iota)*sine(3*beta)+76*et[3,-3]**3 *cosine(iota)*sine(iota)*sine(3*beta)-352/3*et[3,-3]**5 *cosine(iota)*sine(iota)*sine(3*beta)
        Cx05[4,-3]=-(625/16)*et[4,-3]**2 *cosine(iota)*sine(iota)*sine(3*beta)+5625/32*et[4,-3]**4 *cosine(iota)*sine(iota)*sine(3*beta)-(1746875*et[4,-3]**6 *cosine(iota)*sine(iota)*sine(3*beta))/6144
        Cx05[5,-3]=-81*et[5,-3]**3 *cosine(iota)*sine(iota)*sine(3*beta)+729/2*et[5,-3]**5 *cosine(iota)*sine(iota)*sine(3*beta)
        Cx05[6,-3]=-(117649/768)*et[6,-3]**4 *cosine(iota)*sine(iota)*sine(3*beta)+(2705927*et[6,-3]**6 *cosine(iota)*sine(iota)*sine(3*beta))/3840
        Cx05[7,-3]=-(4096/15)*et[7,-3]**5 *cosine(iota)*sine(iota)*sine(3*beta)
        Cx05[8,-3]=-((4782969*et[8,-3]**6 *cosine(iota)*sine(iota)*sine(3*beta))/10240)        
        
        ###################scross05_start####################
        Sx05 = np.zeros((10,5))
        
        Sx05[0,1]=-(1/16)*et[0,1]**2*cosine(iota)*cosine(beta)*sine(iota)-5/96*et[0,1]**4*cosine(iota)*cosine(beta)*sine(iota)-(187*et[0,1]**6*cosine(iota)*cosine(beta)*sine(iota))/6144
        Sx05[1,1]=-(1/24)*et[1,1]**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[2,1]=9/256*et[2,1]**4*cosine(iota)*cosine(beta)*sine(iota)-(63*et[2,1]**6*cosine(iota)*cosine(beta)*sine(iota))/1280
        Sx05[3,1]=1/15*et[3,1]**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[4,1]=(625*et[4,1]**6*cosine(iota)*cosine(beta)*sine(iota))/6144
        
        
        Sx05[0,-1]=-(3/2)*cosine(iota)*cosine(beta)*sine(iota)+et[0,-1]**2*cosine(iota)*cosine(beta)*sine(iota)-5/128*et[0,-1]**4*cosine(iota)*cosine(beta)*sine(iota)+13/384*et[0,-1]**6*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[1,-1]=-2*et[1,-1]*cosine(iota)*cosine(beta)*sine(iota)+2*et[1,-1]**3*cosine(iota)*cosine(beta)*sine(iota)-11/24*et[1,-1]**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[2,-1]=-(45/16)*et[2,-1]**2*cosine(iota)*cosine(beta)*sine(iota)+117/32*et[2,-1]**4*cosine(iota)*cosine(beta)*sine(iota)-(2961*et[2,-1]**6*cosine(iota)*cosine(beta)*sine(iota))/2048
        Sx05[3,-1]=-4*et[3,-1]**3*cosine(iota)*cosine(beta)*sine(iota)+19/3*et[3,-1]**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[4,-1]=-(4375/768)*et[4,-1]**4*cosine(iota)*cosine(beta)*sine(iota)+8125/768*et[4,-1]**6*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[5,-1]=-(81/10)*et[5,-1]**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx05[6,-1]=-((117649*et[6,-1]**6*cosine(iota)*cosine(beta)*sine(iota))/10240)
        
        
        Sx05[0,3]=25/256*et[0,3]**4*cosine(iota)*cosine(3*beta)*sine(iota)+(179*et[0,3]**6*cosine(iota)*cosine(3*beta)*sine(iota))/3840
        Sx05[1,3]=13/120*et[1,3]**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx05[2,3]=(1233*et[2,3]**6*cosine(iota)*cosine(3*beta)*sine(iota))/10240  
        
        
        Sx05[0,-3]=19/16*et[0,-3]**2*cosine(iota)*cosine(3*beta)*sine(iota)-65/96*et[0,-3]**4*cosine(iota)*cosine(3*beta)*sine(iota)-(19*et[0,-3]**6*cosine(iota)*cosine(3*beta)*sine(iota))/6144
        Sx05[1,-3]=-6*et[1,-3]*cosine(iota)*cosine(3*beta)*sine(iota)+11*et[1,-3]**3*cosine(iota)*cosine(3*beta)*sine(iota)-133/24*et[1,-3]**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx05[2,-3]=9/2*cosine(iota)*cosine(3*beta)*sine(iota)-27*et[2,-3]**2*cosine(iota)*cosine(3*beta)*sine(iota)+5319/128*et[2,-3]**4*cosine(iota)*cosine(3*beta)*sine(iota)-3141/128*et[2,-3]**6*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx05[3,-3]=16*et[3,-3]*cosine(iota)*cosine(3*beta)*sine(iota)-76*et[3,-3]**3*cosine(iota)*cosine(3*beta)*sine(iota)+352/3*et[3,-3]**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx05[4,-3]=625/16*et[4,-3]**2*cosine(iota)*cosine(3*beta)*sine(iota)-5625/32*et[4,-3]**4*cosine(iota)*cosine(3*beta)*sine(iota)+(1746875*et[4,-3]**6*cosine(iota)*cosine(3*beta)*sine(iota))/6144
        Sx05[5,-3]=81*et[5,-3]**3*cosine(iota)*cosine(3*beta)*sine(iota)-729/2*et[5,-3]**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx05[6,-3]=117649/768*et[6,-3]**4*cosine(iota)*cosine(3*beta)*sine(iota)-(2705927*et[6,-3]**6*cosine(iota)*cosine(3*beta)*sine(iota))/3840
        Sx05[7,-3]=4096/15*et[7,-3]**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx05[8,-3]=(4782969*et[8,-3]**6*cosine(iota)*cosine(3*beta)*sine(iota))/10240        

        
        ###################cplus1_start####################
        Cp1 = np.zeros((10,7))

        Cp1[0,2]=et[0,2]**5*((11341*cosine(2*beta))/4608-(1561*eta*cosine(2*beta))/4608+1243/512*cosine(iota)**2*cosine(2*beta)-(1099*eta*cosine(iota)**2*cosine(2*beta))/4608-13/576*cosine(2*beta)*sine(iota)**2+13/192*eta*cosine(2*beta)*sine(iota)**2-13/576*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+13/192*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[0,2]**3*(913/288*cosine(2*beta)-269/288*eta*cosine(2*beta)+91/32*cosine(iota)**2*cosine(2*beta)+13/288*eta*cosine(iota)**2*cosine(2*beta)-13/144*cosine(2*beta)*sine(iota)**2+13/48*eta*cosine(2*beta)*sine(iota)**2-13/144*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+13/48*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[1,2]=et[1,2]**6*(137/80*cosine(2*beta)-53/288*eta*cosine(2*beta)+419/240*cosine(iota)**2*cosine(2*beta)-(409*eta*cosine(iota)**2*cosine(2*beta))/1440)+et[1,2]**4*(29/9*cosine(2*beta)-17/16*eta*cosine(2*beta)+17/6*cosine(iota)**2*cosine(2*beta)+5/48*eta*cosine(iota)**2*cosine(2*beta)-5/72*cosine(2*beta)*sine(iota)**2+5/24*eta*cosine(2*beta)*sine(iota)**2-5/72*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+5/24*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[2,2]=et[2,2]**5*(1917/512*cosine(2*beta)-(3297*eta*cosine(2*beta))/2560+(8343*cosine(iota)**2*cosine(2*beta))/2560+(429*eta*cosine(iota)**2*cosine(2*beta))/2560-(81*cosine(2*beta)*sine(iota)**2)/1280+(243*eta*cosine(2*beta)*sine(iota)**2)/1280-(81*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)/1280+(243*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)/1280)
        Cp1[3,2]=et[3,2]**6*((4903*cosine(2*beta))/1080-173/108*eta*cosine(2*beta)+157/40*cosine(iota)**2*cosine(2*beta)+131/540*eta*cosine(iota)**2*cosine(2*beta)-17/270*cosine(2*beta)*sine(iota)**2+17/90*eta*cosine(2*beta)*sine(iota)**2-17/270*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+17/90*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        
        Cp1[0,-2]=et[0,-2]**5*(-((3631*cosine(2*beta))/2304)+(1543*eta*cosine(2*beta))/2304-947/768*cosine(iota)**2*cosine(2*beta)-(827*eta*cosine(iota)**2*cosine(2*beta))/2304+(79*cosine(2*beta)*sine(iota)**2)/1152-79/384*eta*cosine(2*beta)*sine(iota)**2+(79*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)/1152-79/384*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[0,-2]**3*(-(121/32)*cosine(2*beta)-109/96*eta*cosine(2*beta)-179/32*cosine(iota)**2*cosine(2*beta)+413/96*eta*cosine(iota)**2*cosine(2*beta)-1/8*cosine(2*beta)*sine(iota)**2+3/8*eta*cosine(2*beta)*sine(iota)**2-1/8*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+3/8*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[0,-2]*(145/12*cosine(2*beta)+5/4*eta*cosine(2*beta)+53/4*cosine(iota)**2*cosine(2*beta)-9/4*eta*cosine(iota)**2*cosine(2*beta)-1/6*cosine(2*beta)*sine(iota)**2+1/2*eta*cosine(2*beta)*sine(iota)**2-1/6*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+1/2*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[1,-2]= 17/3*cosine(2*beta)-13/3*eta*cosine(2*beta)+3*cosine(iota)**2*cosine(2*beta)+11/3*eta*cosine(iota)**2*cosine(2*beta)+2/3*cosine(2*beta)*sine(iota)**2-2*eta*cosine(2*beta)*sine(iota)**2+2/3*cosine(iota)**2*cosine(2*beta)*sine(iota)**2-2*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+et[1,-2]**2*(235/6*cosine(2*beta)+89/6*eta*cosine(2*beta)+105/2*cosine(iota)**2*cosine(2*beta)-151/6*eta*cosine(iota)**2*cosine(2*beta)+2/3*cosine(2*beta)*sine(iota)**2-2*eta*cosine(2*beta)*sine(iota)**2+2/3*cosine(iota)**2*cosine(2*beta)*sine(iota)**2-2*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[1,-2]**6*(2843/864*cosine(2*beta)+4085/864*eta*cosine(2*beta)+707/96*cosine(iota)**2*cosine(2*beta)-6475/864*eta*cosine(iota)**2*cosine(2*beta)+35/54*cosine(2*beta)*sine(iota)**2-35/18*eta*cosine(2*beta)*sine(iota)**2+35/54*cosine(iota)**2*cosine(2*beta)*sine(iota)**2-35/18*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[1,-2]**4*(-(437/48)*cosine(2*beta)-683/48*eta*cosine(2*beta)-367/16*cosine(iota)**2*cosine(2*beta)+1309/48*eta*cosine(iota)**2*cosine(2*beta)-43/24*cosine(2*beta)*sine(iota)**2+43/8*eta*cosine(2*beta)*sine(iota)**2-43/24*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+43/8*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[2,-2]=et[2,-2]*(-(9/4)*cosine(2*beta)-75/4*eta*cosine(2*beta)-63/4*cosine(iota)**2*cosine(2*beta)+87/4*eta*cosine(iota)**2*cosine(2*beta))+et[2,-2]**3*(2511/32*cosine(2*beta)+1821/32*eta*cosine(2*beta)+4077/32*cosine(iota)**2*cosine(2*beta)-2877/32*eta*cosine(iota)**2*cosine(2*beta)+81/16*cosine(2*beta)*sine(iota)**2-243/16*eta*cosine(2*beta)*sine(iota)**2+81/16*cosine(iota)**2*cosine(2*beta)*sine(iota)**2-243/16*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[2,-2]**5*(-(8379/256)*cosine(2*beta)-14793/256*eta*cosine(2*beta)-21717/256*cosine(iota)**2*cosine(2*beta)+25221/256*eta*cosine(iota)**2*cosine(2*beta)-243/32*cosine(2*beta)*sine(iota)**2+729/32*eta*cosine(2*beta)*sine(iota)**2-243/32*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+729/32*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[3,-2]=et[3,-2]**4*(1294/9*cosine(2*beta)+458/3*eta*cosine(2*beta)+818/3*cosine(iota)**2*cosine(2*beta)-234*eta*cosine(iota)**2*cosine(2*beta)+148/9*cosine(2*beta)*sine(iota)**2-148/3*eta*cosine(2*beta)*sine(iota)**2+148/9*cosine(iota)**2*cosine(2*beta)*sine(iota)**2-148/3*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[3,-2]**2*(-(62/3)*cosine(2*beta)-148/3*eta*cosine(2*beta)-58*cosine(iota)**2*cosine(2*beta)+188/3*eta*cosine(iota)**2*cosine(2*beta)-8/3*cosine(2*beta)*sine(iota)**2+8*eta*cosine(2*beta)*sine(iota)**2-8/3*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+8*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[3,-2]**6*(-(2047/24)*cosine(2*beta)-6161/36*eta*cosine(2*beta)-5639/24*cosine(iota)**2*cosine(2*beta)+10003/36*eta*cosine(iota)**2*cosine(2*beta)-139/6*cosine(2*beta)*sine(iota)**2+139/2*eta*cosine(2*beta)*sine(iota)**2-139/6*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+139/2*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[4,-2]=et[4,-2]**5*((1160825*cosine(2*beta))/4608+(1594375*eta*cosine(2*beta))/4608+277175/512*cosine(iota)**2*cosine(2*beta)-(2406875*eta*cosine(iota)**2*cosine(2*beta))/4608+(94375*cosine(2*beta)*sine(iota)**2)/2304-94375/768*eta*cosine(2*beta)*sine(iota)**2+(94375*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)/2304-94375/768*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[4,-2]**3*(-(15175/288)*cosine(2*beta)-30625/288*eta*cosine(2*beta)-4325/32*cosine(iota)**2*cosine(2*beta)+40625/288*eta*cosine(iota)**2*cosine(2*beta)-625/72*cosine(2*beta)*sine(iota)**2+625/24*eta*cosine(2*beta)*sine(iota)**2-625/72*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+625/24*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[5,-2]=et[5,-2]**6*(68409/160*cosine(2*beta)+11367/16*eta*cosine(2*beta)+163017/160*cosine(iota)**2*cosine(2*beta)-85077/80*eta*cosine(iota)**2*cosine(2*beta)+891/10*cosine(2*beta)*sine(iota)**2-2673/10*eta*cosine(2*beta)*sine(iota)**2+891/10*cosine(iota)**2*cosine(2*beta)*sine(iota)**2-2673/10*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)+et[5,-2]**4*(-(1665/16)*cosine(2*beta)-1647/8*eta*cosine(2*beta)-4257/16*cosine(iota)**2*cosine(2*beta)+2241/8*eta*cosine(iota)**2*cosine(2*beta)-81/4*cosine(2*beta)*sine(iota)**2+243/4*eta*cosine(2*beta)*sine(iota)**2-81/4*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+243/4*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[6,-2]=et[6,-2]**5*(-((844711*cosine(2*beta))/4608)-(8588377*eta*cosine(2*beta))/23040-(3682399*cosine(iota)**2*cosine(2*beta))/7680+(11882549*eta*cosine(iota)**2*cosine(2*beta))/23040-(117649*cosine(2*beta)*sine(iota)**2)/2880+117649/960*eta*cosine(2*beta)*sine(iota)**2-(117649*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)/2880+117649/960*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)
        Cp1[7,-2]=et[7,-2]**6*(-(81713/270)*cosine(2*beta)-17408/27*eta*cosine(2*beta)-24553/30*cosine(iota)**2*cosine(2*beta)+121856/135*eta*cosine(iota)**2*cosine(2*beta)-2048/27*cosine(2*beta)*sine(iota)**2+2048/9*eta*cosine(2*beta)*sine(iota)**2-2048/27*cosine(iota)**2*cosine(2*beta)*sine(iota)**2+2048/9*eta*cosine(iota)**2*cosine(2*beta)*sine(iota)**2)

        
        
        Cp1[0,4]=et[0,4]**5*((1091*cosine(4*beta)*sine(iota)**2)/92160-(1091*eta*cosine(4*beta)*sine(iota)**2)/30720+(1091*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/92160-(1091*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/30720)
        Cp1[1,4]=et[1,4]**6*((71*cosine(4*beta)*sine(iota)**2)/4320-(71*eta*cosine(4*beta)*sine(iota)**2)/1440+(71*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/4320-(71*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/1440)

      
    
        Cp1[0,-4]=et[0,-4]**3*((187*cosine(4*beta)*sine(iota)**2)/1152-187/384*eta*cosine(4*beta)*sine(iota)**2+(187*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/1152-187/384*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[0,-4]**5*(-((1769*cosine(4*beta)*sine(iota)**2)/18432)+(1769*eta*cosine(4*beta)*sine(iota)**2)/6144-(1769*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/18432+(1769*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/6144)
        Cp1[1,-4]=et[1,-4]**4*(137/36*cosine(4*beta)*sine(iota)**2-137/12*eta*cosine(4*beta)*sine(iota)**2+137/36*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-137/12*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[1,-4]**6*(-(57/32)*cosine(4*beta)*sine(iota)**2+171/32*eta*cosine(4*beta)*sine(iota)**2-57/32*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+171/32*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[1,-4]**2*(-(7/3)*cosine(4*beta)*sine(iota)**2+7*eta*cosine(4*beta)*sine(iota)**2-7/3*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+7*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)
        Cp1[2,-4]=et[2,-4]**5*((25191*cosine(4*beta)*sine(iota)**2)/1024-(75573*eta*cosine(4*beta)*sine(iota)**2)/1024+(25191*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/1024-(75573*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/1024)+et[2,-4]*(81/16*cosine(4*beta)*sine(iota)**2-243/16*eta*cosine(4*beta)*sine(iota)**2+81/16*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-243/16*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[2,-4]**3*(-(2511/128)*cosine(4*beta)*sine(iota)**2+7533/128*eta*cosine(4*beta)*sine(iota)**2-2511/128*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+7533/128*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)
        Cp1[3,-4]=-(8/3)*cosine(4*beta)*sine(iota)**2+8*eta*cosine(4*beta)*sine(iota)**2-8/3*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+8*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+et[3,-4]**6*(2768/27*cosine(4*beta)*sine(iota)**2-2768/9*eta*cosine(4*beta)*sine(iota)**2+2768/27*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-2768/9*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[3,-4]**2*(88/3*cosine(4*beta)*sine(iota)**2-88*eta*cosine(4*beta)*sine(iota)**2+88/3*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-88*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[3,-4]**4*(-(253/3)*cosine(4*beta)*sine(iota)**2+253*eta*cosine(4*beta)*sine(iota)**2-253/3*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+253*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)
        Cp1[4,-4]=et[4,-4]**3*(13125/128*cosine(4*beta)*sine(iota)**2-39375/128*eta*cosine(4*beta)*sine(iota)**2+13125/128*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-39375/128*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[4,-4]*(-(625/48)*cosine(4*beta)*sine(iota)**2+625/16*eta*cosine(4*beta)*sine(iota)**2-625/48*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+625/16*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[4,-4]**5*(-((2493125*cosine(4*beta)*sine(iota)**2)/9216)+(2493125*eta*cosine(4*beta)*sine(iota)**2)/3072-(2493125*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/9216+(2493125*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/3072)
        Cp1[5,-4]=et[5,-4]**4*(567/2*cosine(4*beta)*sine(iota)**2-1701/2*eta*cosine(4*beta)*sine(iota)**2+567/2*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-1701/2*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[5,-4]**2*(-(81/2)*cosine(4*beta)*sine(iota)**2+243/2*eta*cosine(4*beta)*sine(iota)**2-81/2*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+243/2*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[5,-4]**6*(-(11745/16)*cosine(4*beta)*sine(iota)**2+35235/16*eta*cosine(4*beta)*sine(iota)**2-11745/16*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+35235/16*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)
        Cp1[6,-4]=et[6,-4]**5*((12588443*cosine(4*beta)*sine(iota)**2)/18432-(12588443*eta*cosine(4*beta)*sine(iota)**2)/6144+(12588443*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/18432-(12588443*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/6144)+et[6,-4]**3*(-((117649*cosine(4*beta)*sine(iota)**2)/1152)+117649/384*eta*cosine(4*beta)*sine(iota)**2-(117649*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/1152+117649/384*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)
        Cp1[7,-4]=et[7,-4]**6*(22528/15*cosine(4*beta)*sine(iota)**2-22528/5*eta*cosine(4*beta)*sine(iota)**2+22528/15*cosine(iota)**2*cosine(4*beta)*sine(iota)**2-22528/5*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)+et[7,-4]**4*(-(2048/9)*cosine(4*beta)*sine(iota)**2+2048/3*eta*cosine(4*beta)*sine(iota)**2-2048/9*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+2048/3*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)
        Cp1[8,-4]=et[8,-4]**5*(-((4782969*cosine(4*beta)*sine(iota)**2)/10240)+(14348907*eta*cosine(4*beta)*sine(iota)**2)/10240-(4782969*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/10240+(14348907*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)/10240)
        Cp1[9,-4]=et[9,-4]**6*(-(390625/432)*cosine(4*beta)*sine(iota)**2+390625/144*eta*cosine(4*beta)*sine(iota)**2-390625/432*cosine(iota)**2*cosine(4*beta)*sine(iota)**2+390625/144*eta*cosine(iota)**2*cosine(4*beta)*sine(iota)**2)    
        
        Cp1[0,0]=et[0,0]**3*(-(11/16)+(23*eta)/48+(11*cosine(iota)**2)/16-23/48*eta*cosine(iota)**2+(9*sine(iota)**2)/64-27/64*eta*sine(iota)**2+9/64*cosine(iota)**2*sine(iota)**2-27/64*eta*cosine(iota)**2*sine(iota)**2)+et[0,0]*(-(9/2)+eta/6+(9*cosine(iota)**2)/2-1/6*eta*cosine(iota)**2+sine(iota)**2/8-3/8*eta*sine(iota)**2+1/8*cosine(iota)**2*sine(iota)**2-3/8*eta*cosine(iota)**2*sine(iota)**2)+et[0,0]**5*(-(235/128)-(47*eta)/1152+(235*cosine(iota)**2)/128+(47*eta*cosine(iota)**2)/1152-(19*sine(iota)**2)/1536+19/512*eta*sine(iota)**2-(19*cosine(iota)**2*sine(iota)**2)/1536+19/512*eta*cosine(iota)**2*sine(iota)**2)
        Cp1[1,0]=et[1,0]**4*(5/2+(35*eta)/18-(5*cosine(iota)**2)/2-35/18*eta*cosine(iota)**2+(7*sine(iota)**2)/12-7/4*eta*sine(iota)**2+7/12*cosine(iota)**2*sine(iota)**2-7/4*eta*cosine(iota)**2*sine(iota)**2)+et[1,0]**6*(-(103/48)-(59*eta)/144+(103*cosine(iota)**2)/48+59/144*eta*cosine(iota)**2-sine(iota)**2/8+3/8*eta*sine(iota)**2-1/8*cosine(iota)**2*sine(iota)**2+3/8*eta*cosine(iota)**2*sine(iota)**2)+et[1,0]**2*(-(15/2)-(11*eta)/6+(15*cosine(iota)**2)/2+11/6*eta*cosine(iota)**2-sine(iota)**2/2+3/2*eta*sine(iota)**2-1/2*cosine(iota)**2*sine(iota)**2+3/2*eta*cosine(iota)**2*sine(iota)**2)
        Cp1[2,0]=et[2,0]**5*(2097/256+(1269*eta)/256-(2097*cosine(iota)**2)/256-1269/256*eta*cosine(iota)**2+(1539*sine(iota)**2)/1024-(4617*eta*sine(iota)**2)/1024+(1539*cosine(iota)**2*sine(iota)**2)/1024-(4617*eta*cosine(iota)**2*sine(iota)**2)/1024)+et[2,0]**3*(-(189/16)-(69*eta)/16+(189*cosine(iota)**2)/16+69/16*eta*cosine(iota)**2-(81*sine(iota)**2)/64+243/64*eta*sine(iota)**2-81/64*cosine(iota)**2*sine(iota)**2+243/64*eta*cosine(iota)**2*sine(iota)**2)
        Cp1[3,0]=et[3,0]**6*(272/15+(472*eta)/45-(272*cosine(iota)**2)/15-472/45*eta*cosine(iota)**2+(16*sine(iota)**2)/5-48/5*eta*sine(iota)**2+16/5*cosine(iota)**2*sine(iota)**2-48/5*eta*cosine(iota)**2*sine(iota)**2)+et[3,0]**4*(-18-(70*eta)/9+18*cosine(iota)**2+70/9*eta*cosine(iota)**2-(7*sine(iota)**2)/3+7*eta*sine(iota)**2-7/3*cosine(iota)**2*sine(iota)**2+7*eta*cosine(iota)**2*sine(iota)**2)
        Cp1[4,0]=et[4,0]**5*(-(6875/256)-(29375*eta)/2304+(6875*cosine(iota)**2)/256+(29375*eta*cosine(iota)**2)/2304-(11875*sine(iota)**2)/3072+(11875*eta*sine(iota)**2)/1024-(11875*cosine(iota)**2*sine(iota)**2)/3072+(11875*eta*cosine(iota)**2*sine(iota)**2)/1024)
        Cp1[5,0]=et[5,0]**6*(-(3159/80)-(1593*eta)/80+(3159*cosine(iota)**2)/80+1593/80*eta*cosine(iota)**2-(243*sine(iota)**2)/40+729/40*eta*sine(iota)**2-243/40*cosine(iota)**2*sine(iota)**2+729/40*eta*cosine(iota)**2*sine(iota)**2)
                    
        ###################splus1_start####################
        Sp1 = np.zeros((10,7))

        Sp1[0,2]=et[0,2]**3*(-(913/288)*sine(2*beta)+269/288*eta*sine(2*beta)-91/32*cosine(iota)**2*sine(2*beta)-13/288*eta*cosine(iota)**2*sine(2*beta)+13/144*sine(iota)**2*sine(2*beta)-13/48*eta*sine(iota)**2*sine(2*beta)+13/144*cosine(iota)**2*sine(iota)**2*sine(2*beta)-13/48*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[0,2]**5*(-((11341*sine(2*beta))/4608)+(1561*eta*sine(2*beta))/4608-1243/512*cosine(iota)**2*sine(2*beta)+(1099*eta*cosine(iota)**2*sine(2*beta))/4608+13/576*sine(iota)**2*sine(2*beta)-13/192*eta*sine(iota)**2*sine(2*beta)+13/576*cosine(iota)**2*sine(iota)**2*sine(2*beta)-13/192*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[1,2]=et[1,2]**6*(-(137/80)*sine(2*beta)+53/288*eta*sine(2*beta)-419/240*cosine(iota)**2*sine(2*beta)+(409*eta*cosine(iota)**2*sine(2*beta))/1440)+et[1,2]**4*(-(29/9)*sine(2*beta)+17/16*eta*sine(2*beta)-17/6*cosine(iota)**2*sine(2*beta)-5/48*eta*cosine(iota)**2*sine(2*beta)+5/72*sine(iota)**2*sine(2*beta)-5/24*eta*sine(iota)**2*sine(2*beta)+5/72*cosine(iota)**2*sine(iota)**2*sine(2*beta)-5/24*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[2,2]=et[2,2]**5*(-(1917/512)*sine(2*beta)+(3297*eta*sine(2*beta))/2560-(8343*cosine(iota)**2*sine(2*beta))/2560-(429*eta*cosine(iota)**2*sine(2*beta))/2560+(81*sine(iota)**2*sine(2*beta))/1280-(243*eta*sine(iota)**2*sine(2*beta))/1280+(81*cosine(iota)**2*sine(iota)**2*sine(2*beta))/1280-(243*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))/1280)
        Sp1[3,2]=et[3,2]**6*(-((4903*sine(2*beta))/1080)+173/108*eta*sine(2*beta)-157/40*cosine(iota)**2*sine(2*beta)-131/540*eta*cosine(iota)**2*sine(2*beta)+17/270*sine(iota)**2*sine(2*beta)-17/90*eta*sine(iota)**2*sine(2*beta)+17/270*cosine(iota)**2*sine(iota)**2*sine(2*beta)-17/90*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        
        
        Sp1[0,-2]=et[0,-2]**5*(-((3631*sine(2*beta))/2304)+(1543*eta*sine(2*beta))/2304-947/768*cosine(iota)**2*sine(2*beta)-(827*eta*cosine(iota)**2*sine(2*beta))/2304+(79*sine(iota)**2*sine(2*beta))/1152-79/384*eta*sine(iota)**2*sine(2*beta)+(79*cosine(iota)**2*sine(iota)**2*sine(2*beta))/1152-79/384*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[0,-2]**3*(-(121/32)*sine(2*beta)-109/96*eta*sine(2*beta)-179/32*cosine(iota)**2*sine(2*beta)+413/96*eta*cosine(iota)**2*sine(2*beta)-1/8*sine(iota)**2*sine(2*beta)+3/8*eta*sine(iota)**2*sine(2*beta)-1/8*cosine(iota)**2*sine(iota)**2*sine(2*beta)+3/8*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[0,-2]*(145/12*sine(2*beta)+5/4*eta*sine(2*beta)+53/4*cosine(iota)**2*sine(2*beta)-9/4*eta*cosine(iota)**2*sine(2*beta)-1/6*sine(iota)**2*sine(2*beta)+1/2*eta*sine(iota)**2*sine(2*beta)-1/6*cosine(iota)**2*sine(iota)**2*sine(2*beta)+1/2*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[1,-2]= 17/3*sine(2*beta)-13/3*eta*sine(2*beta)+3*cosine(iota)**2*sine(2*beta)+11/3*eta*cosine(iota)**2*sine(2*beta)+2/3*sine(iota)**2*sine(2*beta)-2*eta*sine(iota)**2*sine(2*beta)+2/3*cosine(iota)**2*sine(iota)**2*sine(2*beta)-2*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta)+et[1,-2]**2*(235/6*sine(2*beta)+89/6*eta*sine(2*beta)+105/2*cosine(iota)**2*sine(2*beta)-151/6*eta*cosine(iota)**2*sine(2*beta)+2/3*sine(iota)**2*sine(2*beta)-2*eta*sine(iota)**2*sine(2*beta)+2/3*cosine(iota)**2*sine(iota)**2*sine(2*beta)-2*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[1,-2]**6*(2843/864*sine(2*beta)+4085/864*eta*sine(2*beta)+707/96*cosine(iota)**2*sine(2*beta)-6475/864*eta*cosine(iota)**2*sine(2*beta)+35/54*sine(iota)**2*sine(2*beta)-35/18*eta*sine(iota)**2*sine(2*beta)+35/54*cosine(iota)**2*sine(iota)**2*sine(2*beta)-35/18*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[1,-2]**4*(-(437/48)*sine(2*beta)-683/48*eta*sine(2*beta)-367/16*cosine(iota)**2*sine(2*beta)+1309/48*eta*cosine(iota)**2*sine(2*beta)-43/24*sine(iota)**2*sine(2*beta)+43/8*eta*sine(iota)**2*sine(2*beta)-43/24*cosine(iota)**2*sine(iota)**2*sine(2*beta)+43/8*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[2,-2]=et[2,-2]*(-(9/4)*sine(2*beta)-75/4*eta*sine(2*beta)-63/4*cosine(iota)**2*sine(2*beta)+87/4*eta*cosine(iota)**2*sine(2*beta))+et[2,-2]**3*(2511/32*sine(2*beta)+1821/32*eta*sine(2*beta)+4077/32*cosine(iota)**2*sine(2*beta)-2877/32*eta*cosine(iota)**2*sine(2*beta)+81/16*sine(iota)**2*sine(2*beta)-243/16*eta*sine(iota)**2*sine(2*beta)+81/16*cosine(iota)**2*sine(iota)**2*sine(2*beta)-243/16*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[2,-2]**5*(-(8379/256)*sine(2*beta)-14793/256*eta*sine(2*beta)-21717/256*cosine(iota)**2*sine(2*beta)+25221/256*eta*cosine(iota)**2*sine(2*beta)-243/32*sine(iota)**2*sine(2*beta)+729/32*eta*sine(iota)**2*sine(2*beta)-243/32*cosine(iota)**2*sine(iota)**2*sine(2*beta)+729/32*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[3,-2]=et[3,-2]**4*(1294/9*sine(2*beta)+458/3*eta*sine(2*beta)+818/3*cosine(iota)**2*sine(2*beta)-234*eta*cosine(iota)**2*sine(2*beta)+148/9*sine(iota)**2*sine(2*beta)-148/3*eta*sine(iota)**2*sine(2*beta)+148/9*cosine(iota)**2*sine(iota)**2*sine(2*beta)-148/3*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[3,-2]**2*(-(62/3)*sine(2*beta)-148/3*eta*sine(2*beta)-58*cosine(iota)**2*sine(2*beta)+188/3*eta*cosine(iota)**2*sine(2*beta)-8/3*sine(iota)**2*sine(2*beta)+8*eta*sine(iota)**2*sine(2*beta)-8/3*cosine(iota)**2*sine(iota)**2*sine(2*beta)+8*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[3,-2]**6*(-(2047/24)*sine(2*beta)-6161/36*eta*sine(2*beta)-5639/24*cosine(iota)**2*sine(2*beta)+10003/36*eta*cosine(iota)**2*sine(2*beta)-139/6*sine(iota)**2*sine(2*beta)+139/2*eta*sine(iota)**2*sine(2*beta)-139/6*cosine(iota)**2*sine(iota)**2*sine(2*beta)+139/2*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[4,-2]=et[4,-2]**5*((1160825*sine(2*beta))/4608+(1594375*eta*sine(2*beta))/4608+277175/512*cosine(iota)**2*sine(2*beta)-(2406875*eta*cosine(iota)**2*sine(2*beta))/4608+(94375*sine(iota)**2*sine(2*beta))/2304-94375/768*eta*sine(iota)**2*sine(2*beta)+(94375*cosine(iota)**2*sine(iota)**2*sine(2*beta))/2304-94375/768*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[4,-2]**3*(-(15175/288)*sine(2*beta)-30625/288*eta*sine(2*beta)-4325/32*cosine(iota)**2*sine(2*beta)+40625/288*eta*cosine(iota)**2*sine(2*beta)-625/72*sine(iota)**2*sine(2*beta)+625/24*eta*sine(iota)**2*sine(2*beta)-625/72*cosine(iota)**2*sine(iota)**2*sine(2*beta)+625/24*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[5,-2]=et[5,-2]**6*(68409/160*sine(2*beta)+11367/16*eta*sine(2*beta)+163017/160*cosine(iota)**2*sine(2*beta)-85077/80*eta*cosine(iota)**2*sine(2*beta)+891/10*sine(iota)**2*sine(2*beta)-2673/10*eta*sine(iota)**2*sine(2*beta)+891/10*cosine(iota)**2*sine(iota)**2*sine(2*beta)-2673/10*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))+et[5,-2]**4*(-(1665/16)*sine(2*beta)-1647/8*eta*sine(2*beta)-4257/16*cosine(iota)**2*sine(2*beta)+2241/8*eta*cosine(iota)**2*sine(2*beta)-81/4*sine(iota)**2*sine(2*beta)+243/4*eta*sine(iota)**2*sine(2*beta)-81/4*cosine(iota)**2*sine(iota)**2*sine(2*beta)+243/4*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[6,-2]=et[6,-2]**5*(-((844711*sine(2*beta))/4608)-(8588377*eta*sine(2*beta))/23040-(3682399*cosine(iota)**2*sine(2*beta))/7680+(11882549*eta*cosine(iota)**2*sine(2*beta))/23040-(117649*sine(iota)**2*sine(2*beta))/2880+117649/960*eta*sine(iota)**2*sine(2*beta)-(117649*cosine(iota)**2*sine(iota)**2*sine(2*beta))/2880+117649/960*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))
        Sp1[7,-2]=et[7,-2]**6*(-(81713/270)*sine(2*beta)-17408/27*eta*sine(2*beta)-24553/30*cosine(iota)**2*sine(2*beta)+121856/135*eta*cosine(iota)**2*sine(2*beta)-2048/27*sine(iota)**2*sine(2*beta)+2048/9*eta*sine(iota)**2*sine(2*beta)-2048/27*cosine(iota)**2*sine(iota)**2*sine(2*beta)+2048/9*eta*cosine(iota)**2*sine(iota)**2*sine(2*beta))  
        
        
        Sp1[0,4]=et[0,4]**5*(-((1091*sine(iota)**2*sine(4*beta))/92160)+(1091*eta*sine(iota)**2*sine(4*beta))/30720-(1091*cosine(iota)**2*sine(iota)**2*sine(4*beta))/92160+(1091*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/30720)
        Sp1[1,4]=et[1,4]**6*(-((71*sine(iota)**2*sine(4*beta))/4320)+(71*eta*sine(iota)**2*sine(4*beta))/1440-(71*cosine(iota)**2*sine(iota)**2*sine(4*beta))/4320+(71*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/1440) 
        
        
        Sp1[0,-4]=et[0,-4]**3*((187*sine(iota)**2*sine(4*beta))/1152-187/384*eta*sine(iota)**2*sine(4*beta)+(187*cosine(iota)**2*sine(iota)**2*sine(4*beta))/1152-187/384*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[0,-4]**5*(-((1769*sine(iota)**2*sine(4*beta))/18432)+(1769*eta*sine(iota)**2*sine(4*beta))/6144-(1769*cosine(iota)**2*sine(iota)**2*sine(4*beta))/18432+(1769*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/6144)
        Sp1[1,-4]=et[1,-4]**4*(137/36*sine(iota)**2*sine(4*beta)-137/12*eta*sine(iota)**2*sine(4*beta)+137/36*cosine(iota)**2*sine(iota)**2*sine(4*beta)-137/12*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[1,-4]**6*(-(57/32)*sine(iota)**2*sine(4*beta)+171/32*eta*sine(iota)**2*sine(4*beta)-57/32*cosine(iota)**2*sine(iota)**2*sine(4*beta)+171/32*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[1,-4]**2*(-(7/3)*sine(iota)**2*sine(4*beta)+7*eta*sine(iota)**2*sine(4*beta)-7/3*cosine(iota)**2*sine(iota)**2*sine(4*beta)+7*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))
        Sp1[2,-4]=et[2,-4]**5*((25191*sine(iota)**2*sine(4*beta))/1024-(75573*eta*sine(iota)**2*sine(4*beta))/1024+(25191*cosine(iota)**2*sine(iota)**2*sine(4*beta))/1024-(75573*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/1024)+et[2,-4]*(81/16*sine(iota)**2*sine(4*beta)-243/16*eta*sine(iota)**2*sine(4*beta)+81/16*cosine(iota)**2*sine(iota)**2*sine(4*beta)-243/16*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[2,-4]**3*(-(2511/128)*sine(iota)**2*sine(4*beta)+7533/128*eta*sine(iota)**2*sine(4*beta)-2511/128*cosine(iota)**2*sine(iota)**2*sine(4*beta)+7533/128*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))
        Sp1[3,-4]=-(8/3)*sine(iota)**2*sine(4*beta)+8*eta*sine(iota)**2*sine(4*beta)-8/3*cosine(iota)**2*sine(iota)**2*sine(4*beta)+8*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta)+et[3,-4]**6*(2768/27*sine(iota)**2*sine(4*beta)-2768/9*eta*sine(iota)**2*sine(4*beta)+2768/27*cosine(iota)**2*sine(iota)**2*sine(4*beta)-2768/9*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[3,-4]**2*(88/3*sine(iota)**2*sine(4*beta)-88*eta*sine(iota)**2*sine(4*beta)+88/3*cosine(iota)**2*sine(iota)**2*sine(4*beta)-88*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[3,-4]**4*(-(253/3)*sine(iota)**2*sine(4*beta)+253*eta*sine(iota)**2*sine(4*beta)-253/3*cosine(iota)**2*sine(iota)**2*sine(4*beta)+253*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))
        Sp1[4,-4]=et[4,-4]**3*(13125/128*sine(iota)**2*sine(4*beta)-39375/128*eta*sine(iota)**2*sine(4*beta)+13125/128*cosine(iota)**2*sine(iota)**2*sine(4*beta)-39375/128*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[4,-4]*(-(625/48)*sine(iota)**2*sine(4*beta)+625/16*eta*sine(iota)**2*sine(4*beta)-625/48*cosine(iota)**2*sine(iota)**2*sine(4*beta)+625/16*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[4,-4]**5*(-((2493125*sine(iota)**2*sine(4*beta))/9216)+(2493125*eta*sine(iota)**2*sine(4*beta))/3072-(2493125*cosine(iota)**2*sine(iota)**2*sine(4*beta))/9216+(2493125*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/3072)
        Sp1[5,-4]=et[5,-4]**4*(567/2*sine(iota)**2*sine(4*beta)-1701/2*eta*sine(iota)**2*sine(4*beta)+567/2*cosine(iota)**2*sine(iota)**2*sine(4*beta)-1701/2*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[5,-4]**2*(-(81/2)*sine(iota)**2*sine(4*beta)+243/2*eta*sine(iota)**2*sine(4*beta)-81/2*cosine(iota)**2*sine(iota)**2*sine(4*beta)+243/2*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[5,-4]**6*(-(11745/16)*sine(iota)**2*sine(4*beta)+35235/16*eta*sine(iota)**2*sine(4*beta)-11745/16*cosine(iota)**2*sine(iota)**2*sine(4*beta)+35235/16*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))
        Sp1[6,-4]=et[6,-4]**5*((12588443*sine(iota)**2*sine(4*beta))/18432-(12588443*eta*sine(iota)**2*sine(4*beta))/6144+(12588443*cosine(iota)**2*sine(iota)**2*sine(4*beta))/18432-(12588443*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/6144)+et[6,-4]**3*(-((117649*sine(iota)**2*sine(4*beta))/1152)+117649/384*eta*sine(iota)**2*sine(4*beta)-(117649*cosine(iota)**2*sine(iota)**2*sine(4*beta))/1152+117649/384*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))
        Sp1[7,-4]=et[7,-4]**6*(22528/15*sine(iota)**2*sine(4*beta)-22528/5*eta*sine(iota)**2*sine(4*beta)+22528/15*cosine(iota)**2*sine(iota)**2*sine(4*beta)-22528/5*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))+et[7,-4]**4*(-(2048/9)*sine(iota)**2*sine(4*beta)+2048/3*eta*sine(iota)**2*sine(4*beta)-2048/9*cosine(iota)**2*sine(iota)**2*sine(4*beta)+2048/3*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))
        Sp1[8,-4]=et[8,-4]**5*(-((4782969*sine(iota)**2*sine(4*beta))/10240)+(14348907*eta*sine(iota)**2*sine(4*beta))/10240-(4782969*cosine(iota)**2*sine(iota)**2*sine(4*beta))/10240+(14348907*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))/10240)
        Sp1[9,-4]=et[9,-4]**6*(-(390625/432)*sine(iota)**2*sine(4*beta)+390625/144*eta*sine(iota)**2*sine(4*beta)-390625/432*cosine(iota)**2*sine(iota)**2*sine(4*beta)+390625/144*eta*cosine(iota)**2*sine(iota)**2*sine(4*beta))       
                
        ###################ccross1_start####################
        Cx1 = np.zeros((10,7))

        Cx1[0,2]=et[0,2]**5*(-(44/9)*cosine(iota)*sine(2*beta)+(665*eta*cosine(iota)*sine(2*beta))/1152+(131*cosine(iota)*sine(iota)**2*sine(2*beta))/4608-(131*eta*cosine(iota)*sine(iota)**2*sine(2*beta))/1536)+et[0,2]**3*(-(433/72)*cosine(iota)*sine(2*beta)+8/9*eta*cosine(iota)*sine(2*beta)+5/288*cosine(iota)*sine(iota)**2*sine(2*beta)-5/96*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[1,2]=et[1,2]**6*(-(83/24)*cosine(iota)*sine(2*beta)+337/720*eta*cosine(iota)*sine(2*beta)+1/60*cosine(iota)*sine(iota)**2*sine(2*beta)-1/20*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[1,2]**4*(-(109/18)*cosine(iota)*sine(2*beta)+23/24*eta*cosine(iota)*sine(2*beta)-1/18*cosine(iota)*sine(iota)**2*sine(2*beta)+1/6*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[2,2]=et[2,2]**5*(-(2241/320)*cosine(iota)*sine(2*beta)+717/640*eta*cosine(iota)*sine(2*beta)-(297*cosine(iota)*sine(iota)**2*sine(2*beta))/2560+(891*eta*cosine(iota)*sine(iota)**2*sine(2*beta))/2560)
        Cx1[3,2]=et[3,2]**6*(-(4571/540)*cosine(iota)*sine(2*beta)+367/270*eta*cosine(iota)*sine(2*beta)-49/270*cosine(iota)*sine(iota)**2*sine(2*beta)+49/90*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        
        
        Cx1[0,-2]=et[0,-2]*(-(76/3)*cosine(iota)*sine(2*beta)+eta*cosine(iota)*sine(2*beta)+11/12*cosine(iota)*sine(iota)**2*sine(2*beta)-11/4*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[0,-2]**5*(809/288*cosine(iota)*sine(2*beta)-179/576*eta*cosine(iota)*sine(2*beta)+(79*cosine(iota)*sine(iota)**2*sine(2*beta))/2304-79/768*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[0,-2]**3*(75/8*cosine(iota)*sine(2*beta)-19/6*eta*cosine(iota)*sine(2*beta)-21/32*cosine(iota)*sine(iota)**2*sine(2*beta)+63/32*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[1,-2]=-(26/3)*cosine(iota)*sine(2*beta)+2/3*eta*cosine(iota)*sine(2*beta)-8/3*cosine(iota)*sine(iota)**2*sine(2*beta)+8*eta*cosine(iota)*sine(iota)**2*sine(2*beta)+et[1,-2]**2*(-(275/3)*cosine(iota)*sine(2*beta)+31/3*eta*cosine(iota)*sine(2*beta)+16/3*cosine(iota)*sine(iota)**2*sine(2*beta)-16*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[1,-2]**6*(-(4603/432)*cosine(iota)*sine(2*beta)+1195/432*eta*cosine(iota)*sine(2*beta)+20/27*cosine(iota)*sine(iota)**2*sine(2*beta)-20/9*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[1,-2]**4*(769/24*cosine(iota)*sine(2*beta)-313/24*eta*cosine(iota)*sine(2*beta)-10/3*cosine(iota)*sine(iota)**2*sine(2*beta)+10*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[2,-2]=et[2,-2]**3*(-(1647/8)*cosine(iota)*sine(2*beta)+33*eta*cosine(iota)*sine(2*beta)+459/32*cosine(iota)*sine(iota)**2*sine(2*beta)-1377/32*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[2,-2]*(18*cosine(iota)*sine(2*beta)-3*eta*cosine(iota)*sine(2*beta)-27/4*cosine(iota)*sine(iota)**2*sine(2*beta)+81/4*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[2,-2]**5*(1881/16*cosine(iota)*sine(2*beta)-2607/64*eta*cosine(iota)*sine(2*beta)-2781/256*cosine(iota)*sine(iota)**2*sine(2*beta)+8343/256*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[3,-2]=et[3,-2]**4*(-(3748/9)*cosine(iota)*sine(2*beta)+244/3*eta*cosine(iota)*sine(2*beta)+284/9*cosine(iota)*sine(iota)**2*sine(2*beta)-284/3*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[3,-2]**2*(236/3*cosine(iota)*sine(2*beta)-40/3*eta*cosine(iota)*sine(2*beta)-40/3*cosine(iota)*sine(iota)**2*sine(2*beta)+40*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[3,-2]**6*(1281/4*cosine(iota)*sine(2*beta)-1921/18*eta*cosine(iota)*sine(2*beta)-57/2*cosine(iota)*sine(iota)**2*sine(2*beta)+171/2*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[4,-2]=et[4,-2]**5*(-(456925/576)*cosine(iota)*sine(2*beta)+(203125*eta*cosine(iota)*sine(2*beta))/1152+(289375*cosine(iota)*sine(iota)**2*sine(2*beta))/4608-(289375*eta*cosine(iota)*sine(iota)**2*sine(2*beta))/1536)+et[4,-2]**3*(13525/72*cosine(iota)*sine(2*beta)-625/18*eta*cosine(iota)*sine(2*beta)-6875/288*cosine(iota)*sine(iota)**2*sine(2*beta)+6875/96*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[5,-2]=et[5,-2]**6*(-(115713/80)*cosine(iota)*sine(2*beta)+14121/40*eta*cosine(iota)*sine(2*beta)+2349/20*cosine(iota)*sine(iota)**2*sine(2*beta)-7047/20*eta*cosine(iota)*sine(iota)**2*sine(2*beta))+et[5,-2]**4*(2961/8*cosine(iota)*sine(2*beta)-297/4*eta*cosine(iota)*sine(2*beta)-81/2*cosine(iota)*sine(iota)**2*sine(2*beta)+243/2*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        Cx1[6,-2]=et[6,-2]**5*(477211/720*cosine(iota)*sine(2*beta)-(823543*eta*cosine(iota)*sine(2*beta))/5760-(1529437*cosine(iota)*sine(iota)**2*sine(2*beta))/23040+(1529437*eta*cosine(iota)*sine(iota)**2*sine(2*beta))/7680)
        Cx1[7,-2]=et[7,-2]**6*(30269/27*cosine(iota)*sine(2*beta)-34816/135*eta*cosine(iota)*sine(2*beta)-14336/135*cosine(iota)*sine(iota)**2*sine(2*beta)+14336/45*eta*cosine(iota)*sine(iota)**2*sine(2*beta))
        
        
        Cx1[0,4]=et[0,4]**5*(-((1091*cosine(iota)*sine(iota)**2*sine(4*beta))/46080)+(1091*eta*cosine(iota)*sine(iota)**2*sine(4*beta))/15360)
        Cx1[1,4]=et[1,4]**6*(-((71*cosine(iota)*sine(iota)**2*sine(4*beta))/2160)+71/720*eta*cosine(iota)*sine(iota)**2*sine(4*beta))    
        
        
        Cx1[0,-4]=et[0,-4]**5*((1769*cosine(iota)*sine(iota)**2*sine(4*beta))/9216-(1769*eta*cosine(iota)*sine(iota)**2*sine(4*beta))/3072)+et[0,-4]**3*(-(187/576)*cosine(iota)*sine(iota)**2*sine(4*beta)+187/192*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[1,-4]=et[1,-4]**2*(14/3*cosine(iota)*sine(iota)**2*sine(4*beta)-14*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[1,-4]**6*(57/16*cosine(iota)*sine(iota)**2*sine(4*beta)-171/16*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[1,-4]**4*(-(137/18)*cosine(iota)*sine(iota)**2*sine(4*beta)+137/6*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[2,-4]=et[2,-4]**3*(2511/64*cosine(iota)*sine(iota)**2*sine(4*beta)-7533/64*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[2,-4]*(-(81/8)*cosine(iota)*sine(iota)**2*sine(4*beta)+243/8*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[2,-4]**5*(-(25191/512)*cosine(iota)*sine(iota)**2*sine(4*beta)+75573/512*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[3,-4]= 16/3*cosine(iota)*sine(iota)**2*sine(4*beta)-16*eta*cosine(iota)*sine(iota)**2*sine(4*beta)+et[3,-4]**4*(506/3*cosine(iota)*sine(iota)**2*sine(4*beta)-506*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[3,-4]**2*(-(176/3)*cosine(iota)*sine(iota)**2*sine(4*beta)+176*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[3,-4]**6*(-(5536/27)*cosine(iota)*sine(iota)**2*sine(4*beta)+5536/9*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[4,-4]=et[4,-4]**5*((2493125*cosine(iota)*sine(iota)**2*sine(4*beta))/4608-(2493125*eta*cosine(iota)*sine(iota)**2*sine(4*beta))/1536)+et[4,-4]*(625/24*cosine(iota)*sine(iota)**2*sine(4*beta)-625/8*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[4,-4]**3*(-(13125/64)*cosine(iota)*sine(iota)**2*sine(4*beta)+39375/64*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[5,-4]=et[5,-4]**6*(11745/8*cosine(iota)*sine(iota)**2*sine(4*beta)-35235/8*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[5,-4]**2*(81*cosine(iota)*sine(iota)**2*sine(4*beta)-243*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[5,-4]**4*(-567*cosine(iota)*sine(iota)**2*sine(4*beta)+1701*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[6,-4]=et[6,-4]**3*(117649/576*cosine(iota)*sine(iota)**2*sine(4*beta)-117649/192*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[6,-4]**5*(-((12588443*cosine(iota)*sine(iota)**2*sine(4*beta))/9216)+(12588443*eta*cosine(iota)*sine(iota)**2*sine(4*beta))/3072)
        Cx1[7,-4]=et[7,-4]**4*(4096/9*cosine(iota)*sine(iota)**2*sine(4*beta)-4096/3*eta*cosine(iota)*sine(iota)**2*sine(4*beta))+et[7,-4]**6*(-(45056/15)*cosine(iota)*sine(iota)**2*sine(4*beta)+45056/5*eta*cosine(iota)*sine(iota)**2*sine(4*beta))
        Cx1[8,-4]=et[8,-4]**5*((4782969*cosine(iota)*sine(iota)**2*sine(4*beta))/5120-(14348907*eta*cosine(iota)*sine(iota)**2*sine(4*beta))/5120)
        Cx1[9,-4]=et[9,-4]**6*(390625/216*cosine(iota)*sine(iota)**2*sine(4*beta)-390625/72*eta*cosine(iota)*sine(iota)**2*sine(4*beta))       
        
        ###################scross1_start####################
        Sx1 = np.zeros((10,7))
        
        Sx1[0,2]=et[0,2]**5*(-(44/9)*cosine(iota)*cosine(2*beta)+(665*eta*cosine(iota)*cosine(2*beta))/1152+(131*cosine(iota)*cosine(2*beta)*sine(iota)**2)/4608-(131*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)/1536)+et[0,2]**3*(-(433/72)*cosine(iota)*cosine(2*beta)+8/9*eta*cosine(iota)*cosine(2*beta)+5/288*cosine(iota)*cosine(2*beta)*sine(iota)**2-5/96*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[1,2]=et[1,2]**6*(-(83/24)*cosine(iota)*cosine(2*beta)+337/720*eta*cosine(iota)*cosine(2*beta)+1/60*cosine(iota)*cosine(2*beta)*sine(iota)**2-1/20*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[1,2]**4*(-(109/18)*cosine(iota)*cosine(2*beta)+23/24*eta*cosine(iota)*cosine(2*beta)-1/18*cosine(iota)*cosine(2*beta)*sine(iota)**2+1/6*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[2,2]=et[2,2]**5*(-(2241/320)*cosine(iota)*cosine(2*beta)+717/640*eta*cosine(iota)*cosine(2*beta)-(297*cosine(iota)*cosine(2*beta)*sine(iota)**2)/2560+(891*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)/2560)
        Sx1[3,2]=et[3,2]**6*(-(4571/540)*cosine(iota)*cosine(2*beta)+367/270*eta*cosine(iota)*cosine(2*beta)-49/270*cosine(iota)*cosine(2*beta)*sine(iota)**2+49/90*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        
        
        Sx1[0,-2]=et[0,-2]**3*(-(75/8)*cosine(iota)*cosine(2*beta)+19/6*eta*cosine(iota)*cosine(2*beta)+21/32*cosine(iota)*cosine(2*beta)*sine(iota)**2-63/32*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[0,-2]**5*(-(809/288)*cosine(iota)*cosine(2*beta)+179/576*eta*cosine(iota)*cosine(2*beta)-(79*cosine(iota)*cosine(2*beta)*sine(iota)**2)/2304+79/768*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[0,-2]*(76/3*cosine(iota)*cosine(2*beta)-eta*cosine(iota)*cosine(2*beta)-11/12*cosine(iota)*cosine(2*beta)*sine(iota)**2+11/4*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[1,-2]= 26/3*cosine(iota)*cosine(2*beta)-2/3*eta*cosine(iota)*cosine(2*beta)+8/3*cosine(iota)*cosine(2*beta)*sine(iota)**2-8*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2+et[1,-2]**4*(-(769/24)*cosine(iota)*cosine(2*beta)+313/24*eta*cosine(iota)*cosine(2*beta)+10/3*cosine(iota)*cosine(2*beta)*sine(iota)**2-10*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[1,-2]**6*(4603/432*cosine(iota)*cosine(2*beta)-1195/432*eta*cosine(iota)*cosine(2*beta)-20/27*cosine(iota)*cosine(2*beta)*sine(iota)**2+20/9*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[1,-2]**2*(275/3*cosine(iota)*cosine(2*beta)-31/3*eta*cosine(iota)*cosine(2*beta)-16/3*cosine(iota)*cosine(2*beta)*sine(iota)**2+16*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[2,-2]=et[2,-2]**5*(-(1881/16)*cosine(iota)*cosine(2*beta)+2607/64*eta*cosine(iota)*cosine(2*beta)+2781/256*cosine(iota)*cosine(2*beta)*sine(iota)**2-8343/256*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[2,-2]*(-18*cosine(iota)*cosine(2*beta)+3*eta*cosine(iota)*cosine(2*beta)+27/4*cosine(iota)*cosine(2*beta)*sine(iota)**2-81/4*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[2,-2]**3*(1647/8*cosine(iota)*cosine(2*beta)-33*eta*cosine(iota)*cosine(2*beta)-459/32*cosine(iota)*cosine(2*beta)*sine(iota)**2+1377/32*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[3,-2]=et[3,-2]**6*(-(1281/4)*cosine(iota)*cosine(2*beta)+1921/18*eta*cosine(iota)*cosine(2*beta)+57/2*cosine(iota)*cosine(2*beta)*sine(iota)**2-171/2*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[3,-2]**2*(-(236/3)*cosine(iota)*cosine(2*beta)+40/3*eta*cosine(iota)*cosine(2*beta)+40/3*cosine(iota)*cosine(2*beta)*sine(iota)**2-40*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[3,-2]**4*(3748/9*cosine(iota)*cosine(2*beta)-244/3*eta*cosine(iota)*cosine(2*beta)-284/9*cosine(iota)*cosine(2*beta)*sine(iota)**2+284/3*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[4,-2]=et[4,-2]**3*(-(13525/72)*cosine(iota)*cosine(2*beta)+625/18*eta*cosine(iota)*cosine(2*beta)+6875/288*cosine(iota)*cosine(2*beta)*sine(iota)**2-6875/96*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[4,-2]**5*(456925/576*cosine(iota)*cosine(2*beta)-(203125*eta*cosine(iota)*cosine(2*beta))/1152-(289375*cosine(iota)*cosine(2*beta)*sine(iota)**2)/4608+(289375*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)/1536)
        Sx1[5,-2]=et[5,-2]**4*(-(2961/8)*cosine(iota)*cosine(2*beta)+297/4*eta*cosine(iota)*cosine(2*beta)+81/2*cosine(iota)*cosine(2*beta)*sine(iota)**2-243/2*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)+et[5,-2]**6*(115713/80*cosine(iota)*cosine(2*beta)-14121/40*eta*cosine(iota)*cosine(2*beta)-2349/20*cosine(iota)*cosine(2*beta)*sine(iota)**2+7047/20*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        Sx1[6,-2]=et[6,-2]**5*(-(477211/720)*cosine(iota)*cosine(2*beta)+(823543*eta*cosine(iota)*cosine(2*beta))/5760+(1529437*cosine(iota)*cosine(2*beta)*sine(iota)**2)/23040-(1529437*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)/7680)
        Sx1[7,-2]=et[7,-2]**6*(-(30269/27)*cosine(iota)*cosine(2*beta)+34816/135*eta*cosine(iota)*cosine(2*beta)+14336/135*cosine(iota)*cosine(2*beta)*sine(iota)**2-14336/45*eta*cosine(iota)*cosine(2*beta)*sine(iota)**2)
        
        
        Sx1[0,4]=et[0,4]**5*(-((1091*cosine(iota)*cosine(4*beta)*sine(iota)**2)/46080)+(1091*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)/15360)
        Sx1[1,4]=et[1,4]**6*(-((71*cosine(iota)*cosine(4*beta)*sine(iota)**2)/2160)+71/720*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2) 
        
        
        Sx1[0,-4]=et[0,-4]**3*(187/576*cosine(iota)*cosine(4*beta)*sine(iota)**2-187/192*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[0,-4]**5*(-((1769*cosine(iota)*cosine(4*beta)*sine(iota)**2)/9216)+(1769*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)/3072)
        Sx1[1,-4]=et[1,-4]**4*(137/18*cosine(iota)*cosine(4*beta)*sine(iota)**2-137/6*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[1,-4]**6*(-(57/16)*cosine(iota)*cosine(4*beta)*sine(iota)**2+171/16*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[1,-4]**2*(-(14/3)*cosine(iota)*cosine(4*beta)*sine(iota)**2+14*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)
        Sx1[2,-4]=et[2,-4]**5*(25191/512*cosine(iota)*cosine(4*beta)*sine(iota)**2-75573/512*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[2,-4]*(81/8*cosine(iota)*cosine(4*beta)*sine(iota)**2-243/8*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[2,-4]**3*(-(2511/64)*cosine(iota)*cosine(4*beta)*sine(iota)**2+7533/64*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)
        Sx1[3,-4]=-(16/3)*cosine(iota)*cosine(4*beta)*sine(iota)**2+16*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2+et[3,-4]**6*(5536/27*cosine(iota)*cosine(4*beta)*sine(iota)**2-5536/9*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[3,-4]**2*(176/3*cosine(iota)*cosine(4*beta)*sine(iota)**2-176*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[3,-4]**4*(-(506/3)*cosine(iota)*cosine(4*beta)*sine(iota)**2+506*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)
        Sx1[4,-4]=et[4,-4]**3*(13125/64*cosine(iota)*cosine(4*beta)*sine(iota)**2-39375/64*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[4,-4]*(-(625/24)*cosine(iota)*cosine(4*beta)*sine(iota)**2+625/8*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[4,-4]**5*(-((2493125*cosine(iota)*cosine(4*beta)*sine(iota)**2)/4608)+(2493125*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)/1536)
        Sx1[5,-4]=et[5,-4]**4*(567*cosine(iota)*cosine(4*beta)*sine(iota)**2-1701*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[5,-4]**2*(-81*cosine(iota)*cosine(4*beta)*sine(iota)**2+243*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[5,-4]**6*(-(11745/8)*cosine(iota)*cosine(4*beta)*sine(iota)**2+35235/8*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)
        Sx1[6,-4]=et[6,-4]**5*((12588443*cosine(iota)*cosine(4*beta)*sine(iota)**2)/9216-(12588443*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)/3072)+et[6,-4]**3*(-(117649/576)*cosine(iota)*cosine(4*beta)*sine(iota)**2+117649/192*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)
        Sx1[7,-4]=et[7,-4]**6*(45056/15*cosine(iota)*cosine(4*beta)*sine(iota)**2-45056/5*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)+et[7,-4]**4*(-(4096/9)*cosine(iota)*cosine(4*beta)*sine(iota)**2+4096/3*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)
        Sx1[8,-4]=et[8,-4]**5*(-((4782969*cosine(iota)*cosine(4*beta)*sine(iota)**2)/5120)+(14348907*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)/5120)
        Sx1[9,-4]=et[9,-4]**6*(-(390625/216)*cosine(iota)*cosine(4*beta)*sine(iota)**2+390625/72*eta*cosine(iota)*cosine(4*beta)*sine(iota)**2)  
        
        
        Sx1[0,0]=et[0,0]*(1/2*cosine(iota)*sine(iota)**2-3/2*eta*cosine(iota)*sine(iota)**2)+et[0,0]**5*(-(11/384)*cosine(iota)*sine(iota)**2+11/128*eta*cosine(iota)*sine(iota)**2)+et[0,0]**3*(-(5/16)*cosine(iota)*sine(iota)**2+15/16*eta*cosine(iota)*sine(iota)**2)
        Sx1[1,0]=et[1,0]**2*(cosine(iota)*sine(iota)**2-3*eta*cosine(iota)*sine(iota)**2)+et[1,0]**6*(1/12*cosine(iota)*sine(iota)**2-1/4*eta*cosine(iota)*sine(iota)**2)+et[1,0]**4*(-(5/6)*cosine(iota)*sine(iota)**2+5/2*eta*cosine(iota)*sine(iota)**2)
        Sx1[2,0]=et[2,0]**3*(27/16*cosine(iota)*sine(iota)**2-81/16*eta*cosine(iota)*sine(iota)**2)+et[2,0]**5*(-(459/256)*cosine(iota)*sine(iota)**2+1377/256*eta*cosine(iota)*sine(iota)**2)
        Sx1[3,0]=et[3,0]**4*(8/3*cosine(iota)*sine(iota)**2-8*eta*cosine(iota)*sine(iota)**2)+et[3,0]**6*(-(52/15)*cosine(iota)*sine(iota)**2+52/5*eta*cosine(iota)*sine(iota)**2)
        Sx1[4,0]=et[4,0]**5*(3125/768*cosine(iota)*sine(iota)**2-3125/256*eta*cosine(iota)*sine(iota)**2)
        Sx1[5,0]=et[5,0]**6*(243/40*cosine(iota)*sine(iota)**2-729/40*eta*cosine(iota)*sine(iota)**2)

   
        
        #fourier_phase
        
        psi = []
        for n in [0,1,2,3,4,-4,-3,-2,-1]: 
            xpsi = np.abs((G*M*2*np.pi*f)/( C**3*(l - (l + n)*k/(1 + k))))**(2/3)        
            
            psi_n = ( - 2*np.pi*f*tc + (l - (l + n)*k/(1 + k))*phic - 1/(256*xpsi**(5/2)*eta)*3*l * (1 +\
        \
         xpsi**(3/2) * ( - 16*np.pi + et0**2 * ( - ((295945*np.pi)/(35088*chi**(28/9))) + (65561*np.pi)/(4080*chi**(19/9))) + et0**4 * ((1968982405*np.pi)/(35961984*chi**(47/9)) - (6211173025*np.pi)/(102085632*chi**(38/9)) - (3048212305*np.pi)/(64000512*chi**(28/9)) + (217859203*np.pi)/(3720960*chi**(19/9))) + et0**6 * ( - ((28409259125*np.pi)/(79847424*chi**(22/3))) + (30628811474315*np.pi)/(97254162432*chi**(19/3)) + (33366234820475*np.pi)/(65594658816*chi**(47/9)) - (20639727962075*np.pi)/(46551048192*chi**(38/9)) - (126468066221755*np.pi)/(846342770688*chi**(28/9)) + (22156798877*np.pi)/(169675776*chi**(19/9)))) +\
             \
              xpsi * ( - (2585/756) - (25*n)/(3*l) + (55*eta)/9 + et0**6 * (( - (213483902125/1117863936) + (14845156625*eta)/39923712)/chi**7 + (223015085012407/2254879424512 + (13164724715*n)/(71320832*l) - (7378552295*eta)/32530432)/chi**(19/3) + (1968906345873305/5969113952256 - (8999675405695*eta)/16398664704)/chi**(44/9) + ( - (862351154377229525/6373211704344576) - (749497416275*n)/(2742610176*l) + (4063675549105*eta)/13134901248)/chi**(38/9) + ( - (2441897241139735/21246121967616) + (9479155594325*eta)/58368466944)/chi**(25/9) + (116789025584125/3112412663808 + (8448925*n)/(99072*l) - (216909251525*eta)/2585060352)/chi**(19/9)) + et0**4 * ((14796093245/503467776 - (1028884705*eta)/17980992)/chi**(44/9) + ( - (259509826776175/13976341456896) - (225548425*n)/(6014496*l) + (1222893635*eta)/28804608)/chi**(38/9) + ( - (14275935425/416003328) + (209699405*eta)/4000032)/chi**(25/9) + (229668231175/13650932736 + (315685*n)/(8256*l) - (426556895*eta)/11337984)/chi**(19/9)) + et0**2 * (( - (2223905/491232) + (154645*eta)/17544)/chi**(25/9) + (69114725/14968128 + (1805*n)/(172*l) - (128365*eta)/12432)/chi**(19/9))) +\
                  \
                   xpsi**2 * ( - (48825515/508032) - (31805*n)/(252*l) + (22105/504 - (10*n)/l)*eta + (3085*eta**2)/72 + et0**6 * (( - (26945014260125/52819070976) + (17350371000625*eta)/6707183616 - (357715525375*eta**2)/119771136)/chi**(23/3) + (631801735840149031/757639486636032 + (37295665117595*n)/(23963799552*l) + ( - (167443372657006267/47352467914752) - (2593450768855*n)/(855849984*l))*eta + (1453574802115*eta**2)/390365184)/chi**7 + (2231629798479149401848091/2016195625924690968576 + (7565563411268608153*n)/(6689612285411328*l) + (3733620624128821835/13379224570822656 + (7557465734101975*n)/(6126018576384*l))*eta - (173415564792655*eta**2)/148551696384)/chi**(19/3) + (43949506831840859555/63177102070677504 - (1344731894414361455*eta)/376054178992128 + (7946157848161165*eta**2)/2066231752704)/chi**(50/9) + ( - (97831425453010633563475/83514566173731323904) - (85028471563286725*n)/(35939163746304*l) + (132713522808330046295/28679452669550592 + (388656700668275*n)/(98733966336*l))*eta - (2107245064767505*eta**2)/472856444928)/chi**(44/9) + ( - (1193205495087383514152561/777567083994785710080) - (333815911878457549195*n)/(205705577776398336*l) + ( - (139449015255055034899/2314187749984481280) - (2180322429577675295*n)/(1695375641014272*l))*eta + (675785495945689*eta**2)/515614740480)/chi**(38/9) + ( - (387035983120116605285/5846592827536441344) + (1095104635088909345*eta)/1338505683959808 - (185468261986684025*eta**2)/191215097708544)/chi**(31/9) + (214994976177874485475/652560888744640512 + (15553485612355*n)/(20771831808*l) + ( - (14009852235779987735/11652873013297152) - (60376787225*n)/(57065472*l))*eta + (1550053258427425*eta**2)/1488994762752)/chi**(25/9) + (194749034542453075/425320718303232 + (546781247071525*n)/(1048057325568*l) + ( - (122213388638482025/1364570637889536) + (61744532676875*n)/(262014331392*l))*eta - (18060683996675*eta**2)/61262936064)/chi**(19/9)) + et0**4 * ((3654447011975/98224939008 - (4300262795285*eta)/18124839936 + (392328884035*eta**2)/1294631424)/chi**(50/9) + ( - (735191339256903775/7044076094275584) - (638978688025*n)/(3031305984*l) + (55579511401449335/125787073112064 + (44433039725*n)/(108260928*l))*eta - (240910046095*eta**2)/518482944)/chi**(44/9) + ( - (359074780345285439107/1705190973672775680) - (100456187745548465*n)/(451108723193856*l) + ( - (41964795442387913/5074973135930880) - (656130734149165*n)/(3717929037312*l))*eta + (203366083643*eta**2)/1130734080)/chi**(38/9) + ( - (94372278903235/7251965779968) + (126823556396665*eta)/733829870592 - (20940952805*eta**2)/93768192)/chi**(31/9) + (1256913822951125/12777273040896 + (1727660975*n)/(7727616*l) + ( - (1182697961961875/3194318260224) - (25377635*n)/(74304*l))*eta + (34290527545*eta**2)/102041856)/chi**(25/9) + (382978332618985/1865441746944 + (1075257552895*n)/(4596742656*l) + ( - (240335362454795/5984958938112) + (121422004625*n)/(1149185664*l))*eta - (35516739065*eta**2)/268697088)/chi**(19/9)) + et0**2 * ((936702035/1485485568 + (3062285*eta)/260064 - (14251675*eta**2)/631584)/chi**(31/9) + (195802015925/15087873024 + (5113565*n)/(173376*l) + ( - (3656612095/67356576) - (355585*n)/(6192*l))*eta + (25287905*eta**2)/447552)/chi**(25/9) + (115250777195/2045440512 + (323580365*n)/(5040288*l) + ( - (72324815665/6562454976) + (36539875*n)/(1260072*l))*eta - (10688155*eta**2)/294624)/chi**(19/9))) +\
                       \
                        et0**6 * ( - (75356125/(3326976*chi**(19/3))) + 17355248095/(455518464*chi**(38/9)) - 1326481225/(101334144*chi**(19/9))) + et0**4 * (5222765/(998944*chi**(38/9)) - 2608555/(444448*chi**(19/9))) - (2355*et0**2)/(1462*chi**(19/9)) +\
                       \
                        xpsi**(5/2) * ((14453*np.pi)/756 - (32*n*np.pi)/l + et0**2 * (( - ((7063901*np.pi)/520128) + (149064749*np.pi*eta)/2210544)/chi**(34/9) + ((26056251325*np.pi)/1077705216 + (680485*n*np.pi)/(12384*l) - (48393605*np.pi*eta)/895104)/chi**(28/9) + ((185734313*np.pi)/4112640 - (12915517*np.pi*eta)/146880)/chi**(25/9) + ( - ((458370775*np.pi)/6837264) - (4909969*n*np.pi)/(46512*l) + (15803101*np.pi*eta)/229824)/chi**(19/9)) + et0**4 * (((14896370333*np.pi)/61544448 - (351697861441*np.pi*eta)/476969472)/chi**(53/9) + ( - ((7525784976509075*np.pi)/38703714803712) - (85031756225*n*np.pi)/(216521856*l) + (461030900395*np.pi*eta)/1036965888)/chi**(47/9) + ( - ((17596253179825*np.pi)/51451158528) + (1223601085925*np.pi*eta)/1837541376)/chi**(44/9) + ((34901256494241693175*np.pi)/79386134731997184 + (84423313781887*n*np.pi)/(193345546752*l) - (15387742160333*np.pi*eta)/39404703744)/chi**(38/9) + ( - ((2408172473789*np.pi)/6790791168) + (992200223893*np.pi*eta)/1697697792)/chi**(34/9) + ((268377522549925*np.pi)/1965734313984 + (368891935*n*np.pi)/(1188864*l) - (498450665645*np.pi*eta)/1632669696)/chi**(28/9) + ((238457223541*np.pi)/696563712 - (17513506613*np.pi*eta)/33488640)/chi**(25/9) + ( - ((1523166085325*np.pi)/6235584768) - (16315826987*n*np.pi)/(42418944*l) + (52513704623*np.pi*eta)/209599488)/chi**(19/9)) + et0**6 * (( - ((34512939466525*np.pi)/13414367232) + (22598442827675*np.pi*eta)/3353591808)/chi**8 + ((6467437465359803*np.pi)/4162854322176 + (4963101217555*n*np.pi)/(1711699968*l) - (2781714215215*np.pi*eta)/780730368)/chi**(22/3) + ((86771422906734395*np.pi)/32677398577152 - (6033875860440055*np.pi*eta)/1167049949184)/chi**7 + ( - ((36873887275009221134023765*np.pi)/12976515386908092923904) - (272900019722212519495*n*np.pi)/(105681129090269184*l) + (1664283962654437115*np.pi*eta)/623334345080832)/chi**(19/3) + ((616055512637722733*np.pi)/132238832173056 - (292997755491718561*np.pi*eta)/33059708043264)/chi**(53/9) + ( - ((1657908371989673247625*np.pi)/917742485425618944) - (1440942051181375*n*np.pi)/(394935865344*l) + (7812596619965525*np.pi*eta)/1891425779712)/chi**(47/9) + ( - ((2341521777112236925*np.pi)/610004935507968) + (10702863543278075*np.pi*eta)/1675837734912)/chi**(44/9) + ((115976875330365146420525*np.pi)/36200077437790715904 + (280538671697210501*n*np.pi)/(88165569318912*l) - (51133467198786559*np.pi*eta)/17968544907264)/chi**(38/9) + ( - ((279594780479556044255*np.pi)/145760537338970112) + (48634782568328640205*np.pi*eta)/19621610795630592)/chi**(34/9) + ((11134784227004313175*np.pi)/25994870568124416 + (805529084215*n*np.pi)/(827449344*l) - (20680348179051695*np.pi*eta)/21590424059904)/chi**(28/9) + ((203940414046321231*np.pi)/177874509496320 - (158334501890329*np.pi*eta)/97733246976)/chi**(25/9) + ( - ((774548060033375*np.pi)/1421713327104) - (8296791966665*n*np.pi)/(9671519232*l) + (26703843023285*np.pi*eta)/47788683264)/chi**(19/9)) + eta * ( - ((65*np.pi)/9) - 65/9*np.pi*(np.log(f) - np.log(l))) - 1675/756*np.pi*(np.log(f) - np.log(l)) - (160*n*np.pi*(np.log(f) - np.log(l)))/(3*l)) +\
                            \
                             xpsi**3 * (13966988843531/4694215680 + (257982425*n)/(508032*l) - (640*np.pi**2)/3 - (6848*gamma)/21 + ( - (20562265315/3048192) - (2393105*n)/(1512*l) + (23575*np.pi**2)/96 + (1845*n*np.pi**2)/(32*l))*eta + (110255/1728 + (475*n)/(24*l))*eta**2 - (127825*eta**3)/1296 - (13696*np.log(2))/21 - (3424*np.log(xpsi))/21 + et0**2 * (( - (82471214720975/45625728024576) - (2153818055*n)/(524289024*l) + ( - (48415393035455/1629490286592) - (119702185*n)/(1560384*l))*eta + (906325428545/6466231296 + (32769775*n)/(222912*l))*eta**2 - (2330466575*eta**3)/16111872)/chi**(31/9) + (24716497*np.pi**2)/(293760*chi**(28/9)) + (326505451793435/2061804036096 + (916703174045*n)/(5080610304*l) + ( - (13467050491570355/39689727694848) - (9519440485*n)/(35282016*l))*eta + ( - (2186530635995/52499639808) - (7198355375*n)/(45362592*l))*eta**2 + (2105566535*eta**3)/10606464)/chi**(25/9) + 1/chi**(19/9)*(4175723876720788380517/5556561877278720000 + (534109712725265*n)/(2405438042112*l) - (21508213*np.pi**2)/276480 - (734341*gamma)/16800 + ( - (37399145056383727/28865256505344) - (1219797059185*n)/(2045440512*l) + (12111605*np.pi**2)/264192 + (639805*n*np.pi**2)/(22016*l))*eta + ( - (159596464273381/1718170030080) + (43766986495*n)/(1022720256*l))*eta**2 - (69237581*eta**3)/746496 - (9663919*np.log(2))/50400 + (4602177*np.log(3))/44800 - (734341*np.log(xpsi))/33600) + 1/chi**(37/9)*( - (4165508390854487/16471063977984) - (96423905*np.pi**2)/5052672 + (2603845*gamma)/61404 + ( - (1437364085977/53477480448) + (3121945*np.pi**2)/561408)*eta + (4499991305*eta**2)/636636672 + (2425890995*eta**3)/68211072 + (1898287*np.log(2))/184212 + (12246471*np.log(3))/163744 + (2603845*np.log(xpsi))/122808 - (2603845*np.log(chi))/184212)) + et0**4 * (1/chi**(50/9)*( - (181582918442691290125/1374276523167055872) - (157819616198875*n)/(591398019072*l) + (1741702918744309017425/1521520436363526144 + (185709581143825*n)/(109127015424*l))*eta + ( - (18130335399490218365/6037779509379072) - (16942972137575*n)/(7794786816*l))*eta**2 + (91862546967565*eta**3)/37330771968) - (2341612230425*np.pi**2)/(3675082752*chi**(47/9)) + 1/chi**(44/9)*( - (1017258852718193648990131/859416250731078942720) - (284592379883138801345*n)/(227358796489703424*l) + (69311096542161812013731/30693437526109962240 + (17602484074819772515*n)/(12179935526234112*l))*eta + (3272123415010135297/2970715982008320 + (129257754627385505*n)/(66922722671616*l))*eta**2 - (40063118477671*eta**3)/20353213440) + 1/chi**(31/9)*(141251897794072110575/3786570420215611392 + (194154433667165*n)/(2290094456832*l) + ( - (11182467092862313645/19319236837834752) - (15348073704055*n)/(13631514624*l))*eta + (1038816664853665/594291769344 + (2534255435*n)/(1741824*l))*eta**2 - (147245442666235*eta**3)/102858190848) + (254578148953*np.pi**2)/(535818240*chi**(28/9)) + 1/chi**(25/9)*(2095939685244436475/1746053475139584 + (5884601777755325*n)/(4302551126016*l) + ( - (17381974915387486205/8402882349109248) - (527634379756765*n)/(358545927168*l))*eta + ( - (386694251193132845/933653594345472) - (9761006428375*n)/(10342670976*l))*eta**2 + (2855158909615*eta**3)/2418273792) + 1/chi**(19/9)*(13875930442343179788457991/5067584432078192640000 + (1774846575386055595*n)/(2193759494406144*l) - (71471791799*np.pi**2)/252149760 - (2440215143*gamma)/15321600 + ( - (124277359022363124821/26325113932873728) - (4053385627671755*n)/(1865441746944*l) + (40246863415*np.pi**2)/240943104 + (2126072015*n*np.pi**2)/(20078592*l))*eta + ( - (530339050780445063/1566971067432960) + (7654615585415*n)/(49090572288*l))*eta**2 - (230076481663*eta**3)/680804352 - (32113202837*np.log(2))/45964800 + (5097678057*np.log(3))/13619200 - (2440215143*np.log(xpsi))/30643200) + 1/chi**(38/9)*( - (3123488330286080905561719773/355085641155718958284800) - (85280660877506238107*n)/(124770071244349440*l) + (300051120571*np.pi**2)/970776576 + (211649317*gamma)/191520 + ( - (40336854286157147692937/32939298808508252160) + (584462420500316711*n)/(495119330334720*l) + (2786391039419*np.pi**2)/17972849664 - (91683875075*n*np.pi**2)/(1089263616*l))*eta + (14654969487690651143/35648591784099840 - (46042929781519*n)/(107385626880*l))*eta**2 + (49171400252465*eta**3)/91738386432 + (2117998887803*np.log(2))/44241120 - (334711679031*np.log(3))/13108480 + (211649317*np.log(xpsi))/383040) + 1/chi**(56/9)*(259620437372696563/159257838845952 + (691917129965*np.pi**2)/2589262848 - (558835855*gamma)/2030112 + ( - (245999063921173/13702378991616) - (20770936405*np.pi**2)/575391744)*eta + (255806950720535*eta**2)/326247118848 - (9022269087085*eta**3)/8738762112 - (12629690323*np.log(2))/188800416 - (27159422553*np.log(3))/55940864 - (558835855*np.log(xpsi))/4060224 + (558835855*np.log(chi))/6090336) + 1/chi**(37/9)*(102453749612934666311/19868699733442560 - (598067688595*np.pi**2)/4608036864 - (36290762107*gamma)/56000448 + (6738669506224179365/2219101528670208 - (110934582115*np.pi**2)/512004096)*eta - (1484623162301215*eta**2)/6604468835328 + (128895671353745*eta**3)/217729741824 - (1140350944327*np.log(2))/24000192 + (1296725746149*np.log(3))/49778176 - (36290762107*np.log(xpsi))/112000896 + (36290762107*np.log(chi))/168001344)) + et0**6 * (1/chi**(23/3)*(79743280932801358583/35798465743552512 + (4707297451617835*n)/(1132289528832*l) + ( - (4697965091339819013485/286387725948420096) - (3031112042005175*n)/(143782797312*l))*eta + (44260937463883607881/1136459229954048 + (62492948222105*n)/(2567549952*l))*eta**2 - (35025987744365*eta**3)/1171095552) + (398174549166095*np.pi**2)/(80486203392*chi**(22/3)) + 1/chi**7 * (6322207219091430255435641803/677441730310696165441536 + (21433241144123966897449*n)/(2247709727898206208*l) + ( - (382703787518451715108530917/24194347511096291622912) - (163853568864047723579*n)/(20068836856233984*l))*eta + ( - (106526264171461811599015/7385331963094106112) - (1488820749618089075*n)/(73512222916608*l))*eta**2 + (34162866264153035*eta**3)/1782620356608) + 1/chi**(50/9)*( - (2183772179071687351051132225/883918168382772332199936) - (1897987377461640667975*n)/(380380109090881536*l) + (579230324539710961384384955/31568506013670440435712 + (58073101284937293475*n)/(2264167316017152*l))*eta + ( - (10549226416696214484085/247084515306897408) - (343159875555216425*n)/(12440479758336*l))*eta**2 + (1860567315439539235*eta**3)/59579912060928) - (39680793155110375*np.pi**2)/(6703350939648*chi**(47/9)) + 1/chi**(44/9)*( - (135365962984146442151489284759/10189239068667671944888320) - (37870519836609632882704205*n)/(2695565891181923794944*l) + (7853835483432029621933284961/363901395309559712317440 + (33994428248774224759285*n)/(2777025299981377536*l))*eta + (1356831692547764227006129/111081011999255101440 + (1130620204248677608895*n)/(61033523076513792*l))*eta**2 - (350432910788522809*eta**3)/18562130657280) + 1/chi**(31/9)*(34076287310129300818097225/179574315608305154654208 + (2465197344714118505*n)/(5716075764252672*l) + ( - (2537511568718810179185985/916195487797475278848) - (6975188758528085*n)/(1308625403904*l))*eta + (2244292522811400055685/279668952319131648 + (1181326509469325*n)/(186946486272*l))*eta**2 - (30328195477605980725*eta**3)/4877946842775552) + (10562258458043923*np.pi**2)/(7085660405760*chi**(28/9)) + 1/chi**(25/9)*(358510260983032848848845/89174443082328834048 + (1006560510293615881915*n)/(219739891107889152*l) + ( - (11101485749253028049541715/1716608029334830055424) - (40317372635313267425*n)/(9155828796162048*l))*eta + ( - (914845583455980713785/619266965849505792) - (441232051620619375*n)/(150920254881792*l))*eta**2 + (129063292052563975*eta**3)/35287451172864) + 1/chi**(19/3)*(580096325546747517323966372962837/10563010179869529959647150080 + (151645106281290420888707167*n)/(45364401583305546792960*l) - (1918746213416491*np.pi**2)/1000328527872 - (5813865129161*gamma)/815109120 + (276968936971755062954790967/16033421036455817379840 - (2022450572459665609*n)/(746979039313920*l) - (7714456204997411*np.pi**2)/5960450506752 + (27993947842265*n*np.pi**2)/(70957744128*l))*eta + ( - (662140667786869733841389/387729928062440570880) + (4063733530370815843*n)/(964041870704640*l))*eta**2 - (3388909956719855*eta**3)/813804945408 - (22635681300089561*np.log(2))/66023838720 + (18242991444087*np.log(3))/103505920 + (910126953125*np.log(5))/146313216 - (5813865129161*np.log(xpsi))/1630218240) + 1/chi**(19/9)*(1411215114204927478714383769/231081850102765584384000 + (902530580917461918025*n)/(500177164724600832*l) - (7268851140841*np.pi**2)/11498029056 - (248175681337*gamma)/698664960 + ( - (63196514329101376128695/6002125976695209984) - (2061194773654925225*n)/(425320718303232*l) + (20466008454925*np.pi**2)/54935027712 + (1081132891925*n*np.pi**2)/(4577918976*l))*eta + ( - (53936742276439022617/71453880674942976) + (3892463014444925*n)/(11192650481664*l))*eta**2 - (116996625810085*eta**3)/155223392256 - (3265989073483*np.log(2))/2095994880 + (172815325821*np.log(3))/207011840 - (248175681337*np.log(xpsi))/1397329920) + 1/chi**(38/9)*( - (10379351721540646849181594805679/161919052367007844977868800) - (283387636095953229229561*n)/(56895152487423344640*l) + (997069873657433*np.pi**2)/442674118656 + (703310680391*gamma)/87333120 + ( - (134039366792900201783629651/15020320256679762984960) + (1942168623322552430653*n)/(225774414632632320*l) + (9259177423989337*np.pi**2)/8195619446784 - (304665516874225*n*np.pi**2)/(496704208896*l))*eta + (48698463607596033748189/16255757853549527040 - (153000655663987637*n)/(48967845857280*l))*eta**2 + (163396563038941195*eta**3)/41832704212992 + (7038110304169369*np.log(2))/20173950720 - (370748969806671*np.log(3))/1992488960 + (703310680391*np.log(xpsi))/174666240) + 1/chi**(56/9)*( - (95765636723679036324982133/3502538538798360821760) + (13756565834952955*np.pi**2)/4722815434752 + (45970619802497*gamma)/14348831616 + ( - (55011254544918787424693/2274375674544390144) + (20907767625235*np.pi**2)/16398664704)*eta + (463380478491174152645*eta**2)/27075900887433216 -\
                             (3879443939044136875*eta**3)/223153029292032 + (13266735591208763*np.log(2))/43046494848 - (726469287588495*np.log(3))/4251505664 + (45970619802497*np.log(xpsi))/28697663232 - (45970619802497*np.log(chi))/43046494848) + 1/chi**(25/3)*( - (1844247076182880525/167330816851968) - (10225600094125*np.pi**2)/3832676352 + (249956266625*gamma)/139732992 + (2728185267249325/633828851712 + (299691309125*np.pi**2)/1277558784)*eta - (1304478350387875*eta**2)/80486203392 + (134782341955625*eta**3)/8623521792 + (182226181475*np.log(2))/419198976 + (130622307075*np.log(3))/41402368 + (249956266625*np.log(xpsi))/279465984 - (249956266625*np.log(chi))/419198976) + 1/chi**(37/9)*(1720927919854684009084595897/40516888294827538513920 - (3986831179520597*np.pi**2)/8405059239936 - (1423526912698421*gamma)/255362042880 + (131543151853096063653535/6190510052685643776 - (775201866281389*np.pi**2)/466947735552)*eta - (174156473319672237061*eta**2)/96372409245106176 + (18782995537481836405*eta**3)/5162807638130688 - (233279096767651103*np.log(2))/766086128640 + (773344207011339*np.log(3))/4450754560 - (164052734375*np.log(5))/20959232 - (1423526912698421*np.log(xpsi))/510724085760 + (1423526912698421*np.log(chi))/766086128640)))))
        
            psi.append(psi_n)
        
        # psi( n , l-1 )
        psi = np.array(psi).astype('float64')     
                
        
        # unit( l-1 , n ) 
        unit = np.zeros((10,9))
        for n in [0,1,2,3,4,-4,-3,-2,-1]: 
            for ll in l:
                unit[ll-1,n] = self.unitstep( (ll-(ll+n)*(k[ll-1]/(1+k[ll-1]))) ,ff,f)   

        ####################################################################################
        #for plus polarization
        xii = Xi(1.0, 0.0, eta, et, unit, Cp0, Sp0, Cx0, Sx0, Cp05, Sp05, Cx05, Sx05, Cp1, Sp1, Cx1, Sx1)
        
        ##########0PN#############
        # xi0( l-1 , n )
        xi0 = xii.xi0_
        
        s1 = 0
        s2 = 0
        s3 = 0                                                                                   
        #frequency domain waveform 
        # psi( n , l-1 )
        for ll in [1,2,3,4,5,6,7,8]:
            #n = -2                                                                                
            s1 = s1 + xi0(ll-1,-2)*((ll/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,-2]) )
            
        for ll in [1,2,3,4]:
            #n = 2                                                                                
            s2 = s2 + xi0(ll-1,2)*((ll/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,2]) )
            
        for ll in [1,2,3,4,5,6]:
            #n = 0                                                                                
            s3 = s3 + xi0(ll-1,0)*((ll/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,0]) ) 
        
        hf0 = s1+s2+s3
        
        ##########05PN#############
        # xi05( l-1 , n )
        xi05 = xii.xi05_
        
        s1 = 0
        s2 = 0
        s3 = 0 
        s4 = 0                                                                                   
        #frequency domain waveform 
        for ll in [1,2,3,4,5,6,7]:
            #n = -1                                                                                
            s1 = s1 + xi05(ll-1,-1)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,-1]) )
        for l in [1,2,3,4,5]:
            #n = 1                                                                                
            s2 = s2 + xi05(ll-1,1)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,1]) )
        for ll in [1,2,3,4,5,6,7,8,9]:
            #n = -3                                                                                
            s3 = s3 + xi05(ll-1,-3)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,-3]) )
        for ll in [1,2,3]:
            #n = 3                                                                                
            s4 = s4 + xi05(ll-1,3)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,3]) )                  
        hf05 = s1+s2+s3+s4
        
        ##########1PN#############
        # xi1( l-1 , n )
        xi1 = xii.xi1_
        xipn = xii.xipn_
        
        s1 = 0
        s2 = 0
        s3 = 0 
        s4 = 0   
        s5 = 0
        #frequency domain waveform 
        for ll in [1,2,3,4,5,6,7,8]:
            #n = -2                                                                                
            s1 = s1 + (xi1(ll-1,-2) + xipn(ll-1,-2))*np.exp( -1j*(np.pi/4 + psi[ll-1,-2]) )
        for ll in [1,2,3,4]:
            #n = 2                                                                                
            s2 = s2 + (xi1(ll-1,2) + xipn(ll-1,2))*np.exp( -1j*(np.pi/4 + psi[ll-1,2]) )
        for ll in [1,2,3,4,5,6]:
            #n = 0                                                                                
            s3 = s3 + (xi1(ll-1,0) + xipn(ll-1,0))*np.exp( -1j*(np.pi/4 + psi[ll-1,0]) )
        for ll in [1,2,3,4,5,6,7,8,9,10]:
            #n = -4                                                                                
            s4 = s4 + xi1(ll-1,-4)*np.exp( -1j*(np.pi/4 + psi[ll-1,-4]) )
        for ll in [1,2]:
            #n = 4                                                                                
            s5 = s5 + xi1(ll-1,4)*np.exp( -1j*(np.pi/4 + psi[ll-1,4]) )
        
        hf1 = s1+s2+s3+s4+s5
 
        # final frequency domain strain
        hp = (((5*np.pi*eta)/384)**(1/2)) * (G**2*M**2)/(C**5*D)*( (((G*M*np.pi*f)/C**3)**(-7/6))*hf0 + (((G*M*np.pi*f)/C**3)**(-5/6))*(delta/M)*hf05 + (((G*M*np.pi*f)/C**3)**(-1/2))*hf1 )
        
        
        
        ####################################################################################
        #for cross polarization
        xii = Xi(0.0, 1.0, eta, et, unit, Cp0, Sp0, Cx0, Sx0, Cp05, Sp05, Cx05, Sx05, Cp1, Sp1, Cx1, Sx1)
        
        ##########0PN#############
        # xi0( l-1 , n )
        xi0 = xii.xi0_
        
        s1 = 0
        s2 = 0
        s3 = 0                                                                                   
        #frequency domain waveform 
        # psi( n , l-1 )
        for ll in [1,2,3,4,5,6,7,8]:
            #n = -2                                                                                
            s1 = s1 + xi0(ll-1,-2)*((ll/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,-2]) )
            
        for ll in [1,2,3,4]:
            #n = 2                                                                                
            s2 = s2 + xi0(ll-1,2)*((ll/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,2]) )
            
        for ll in [1,2,3,4,5,6]:
            #n = 0                                                                                
            s3 = s3 + xi0(ll-1,0)*((ll/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,0]) ) 
        
        hf0 = s1+s2+s3
        
        ##########05PN#############
        # xi05( l-1 , n )
        xi05 = xii.xi05_
        
        s1 = 0
        s2 = 0
        s3 = 0 
        s4 = 0                                                                                   
        #frequency domain waveform 
        for ll in [1,2,3,4,5,6,7]:
            #n = -1                                                                                
            s1 = s1 + xi05(ll-1,-1)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,-1]) )
        for l in [1,2,3,4,5]:
            #n = 1                                                                                
            s2 = s2 + xi05(ll-1,1)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,1]) )
        for ll in [1,2,3,4,5,6,7,8,9]:
            #n = -3                                                                                
            s3 = s3 + xi05(ll-1,-3)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,-3]) )
        for ll in [1,2,3]:
            #n = 3                                                                                
            s4 = s4 + xi05(ll-1,3)*((ll/2)**(1/3))*np.exp( -1j*(np.pi/4 + psi[ll-1,3]) )                  
        hf05 = s1+s2+s3+s4
        
        ##########1PN#############
        # xi1( l-1 , n )
        xi1 = xii.xi1_
        xipn = xii.xipn_
        
        s1 = 0
        s2 = 0
        s3 = 0 
        s4 = 0   
        s5 = 0
        #frequency domain waveform 
        for ll in [1,2,3,4,5,6,7,8]:
            #n = -2                                                                                
            s1 = s1 + (xi1(ll-1,-2) + xipn(ll-1,-2))*np.exp( -1j*(np.pi/4 + psi[ll-1,-2]) )
        for ll in [1,2,3,4]:
            #n = 2                                                                                
            s2 = s2 + (xi1(ll-1,2) + xipn(ll-1,2))*np.exp( -1j*(np.pi/4 + psi[ll-1,2]) )
        for ll in [1,2,3,4,5,6]:
            #n = 0                                                                                
            s3 = s3 + (xi1(ll-1,0) + xipn(ll-1,0))*np.exp( -1j*(np.pi/4 + psi[ll-1,0]) )
        for ll in [1,2,3,4,5,6,7,8,9,10]:
            #n = -4                                                                                
            s4 = s4 + xi1(ll-1,-4)*np.exp( -1j*(np.pi/4 + psi[ll-1,-4]) )
        for ll in [1,2]:
            #n = 4                                                                                
            s5 = s5 + xi1(ll-1,4)*np.exp( -1j*(np.pi/4 + psi[ll-1,4]) )
        
        hf1 = s1+s2+s3+s4+s5
 
        # final frequency domain strain
        hc = (((5*np.pi*eta)/384)**(1/2)) * (G**2*M**2)/(C**5*D)*( (((G*M*np.pi*f)/C**3)**(-7/6))*hf0 + (((G*M*np.pi*f)/C**3)**(-5/6))*(delta/M)*hf05 + (((G*M*np.pi*f)/C**3)**(-1/2))*hf1 )

        return(et)
    
    

class Xi:
    def __init__(self, Fp_, Fc_, eta_, et_, unit_, Cp0_, Sp0_, Cx0_, Sx0_, Cp05_, Sp05_, Cx05_, Sx05_, Cp1_, Sp1_, Cx1_, Sx1_):
        
        self.Fp_ = Fp_
        self.Fc_ = Fc_
        self.eta_ = eta_
        self.et_ = et_
        self.unit_ = unit_
        self.Cp0_ = Cp0_
        self.Cx0_ = Cx0_
        self.Sp0_ = Sp0_
        self.Sx0_ = Sx0_
        self.Cp05_ = Cp05_
        self.Cx05_ = Cx05_
        self.Sp05_ = Sp05_
        self.Sx05_ = Sx05_
        self.Cp1_ = Cp1_
        self.Cx1_ = Cx1_
        self.Sp1_ = Sp1_
        self.Sx1_ = Sx1_
    
    #here l range from 0 to 9
    def xi0_(self, l_, n_):
        Fp = self.Fp_
        Fc = self.Fc_
        et = self.et_
        unit = self.unit_
        Cp0 = self.Cp0_
        Cx0 = self.Cx0_
        Sp0 = self.Sp0_
        Sx0 = self.Sx0_
        
        if unit[l_,n_]==0:
            xil = 0.0*1j
            
        else:    
            Gamma_l = Fp*Cp0[l_,n_] + Fc*Cx0[l_,n_]
            Sigma_l = Fp*Sp0[l_,n_] + Fc*Sx0[l_,n_]

            al = np.sign(Gamma_l)*np.sqrt(Gamma_l**2 + Sigma_l**2)

            if Gamma_l==0:
                phil = -np.sign(Sigma_l)*np.pi/2 
            else:
                phil = np.arctan(- (Sigma_l/Gamma_l))            

            #NOTE: et(n value, l-1)
            numerator = (1-et[l_,n_]**2)**(7/4)
            denomitor = ( 1 + (73/24)*et[l_,n_]**2 + (37/96)*et[l_,n_]**4 )**(1/2)
            xil = (numerator/denomitor)*al*np.exp(-1j*phil)
            
        return(xil)
    
    def xi05_(self, l_, n_):
        Fp = self.Fp_
        Fc = self.Fc_
        et = self.et_
        unit = self.unit_
        Cp05 = self.Cp05_
        Cx05 = self.Cx05_
        Sp05 = self.Sp05_
        Sx05 = self.Sx05_
        
        if unit[l_,n_]==0:
            xil = 0.0*1j
            
        else:    
            Gamma_l = Fp*Cp05[l_,n_] + Fc*Cx05[l_,n_]
            Sigma_l = Fp*Sp05[l_,n_] + Fc*Sx05[l_,n_]

            al = np.sign(Gamma_l)*np.sqrt(Gamma_l**2 + Sigma_l**2)

            if Gamma_l==0:
                phil = -np.sign(Sigma_l)*np.pi/2 
            else:
                phil = np.arctan(- (Sigma_l/Gamma_l))            

            #NOTE: et(n value, l-1)
            numerator = (1-et[l_,n_]**2)**(7/4)
            denomitor = ( 1 + (73/24)*et[l_,n_]**2 + (37/96)*et[l_,n_]**4 )**(1/2)
            xil = (numerator/denomitor)*al*np.exp(-1j*phil)
            
        return(xil)
    
    def xi1_(self, l_, n_):
        Fp = self.Fp_
        Fc = self.Fc_
        et = self.et_
        unit = self.unit_
        Cp1 = self.Cp1_
        Cx1 = self.Cx1_
        Sp1 = self.Sp1_
        Sx1 = self.Sx1_
        
        if unit[l_,n_]==0:
            xil = 0.0*1j
            
        else:    
            Gamma_l = Fp*Cp1[l_,n_] + Fc*Cx1[l_,n_]
            Sigma_l = Fp*Sp1[l_,n_] + Fc*Sx1[l_,n_]

            al = np.sign(Gamma_l)*np.sqrt(Gamma_l**2 + Sigma_l**2)

            if Gamma_l==0:
                phil = -np.sign(Sigma_l)*np.pi/2 
            else:
                phil = np.arctan(- (Sigma_l/Gamma_l))            

            #NOTE: et(n value, l-1)
            numerator = (1-et[l_,n_]**2)**(7/4)
            denomitor = ( 1 + (73/24)*et[l_,n_]**2 + (37/96)*et[l_,n_]**4 )**(1/2)
            xil = (numerator/denomitor)*al*np.exp(-1j*phil)
            
        return(xil)
    
    def xipn_(self, l_, n_):
        Fp = self.Fp_
        Fc = self.Fc_
        eta = self.eta_
        et = self.et_
        unit = self.unit_
        Cp0 = self.Cp0_
        Cx0 = self.Cx0_
        Sp0 = self.Sp0_
        Sx0 = self.Sx0_
        
        if unit[l_,n_]==0:
            xil = 0.0*1j
            
        else:    
            Gamma_l = Fp*Cp0[l_,n_] + Fc*Cx0[l_,n_]
            Sigma_l = Fp*Sp0[l_,n_] + Fc*Sx0[l_,n_]

            al = np.sign(Gamma_l)*np.sqrt(Gamma_l**2 + Sigma_l**2)

            if Gamma_l==0:
                phil = -np.sign(Sigma_l)*np.pi/2 
            else:
                phil = np.arctan(- (Sigma_l/Gamma_l))            

            #NOTE: et(n value, l-1)
            numerator = (1-et[l_,n_]**2)**(3/4)*(11888 + 14784*eta - et[l_,n_]**2*(87720-159600*eta) - et[l_,n_]**4*(171038-141708*eta) - et[l_,n_]**6*(11717-8288*eta))
            denomitor = ( 1 + (73/24)*et[l_,n_]**2 + (37/96)*et[l_,n_]**4 )**(3/2) * 10752

            xil = (numerator/denomitor)*al*np.exp(-1j*phil) 
            
        return(xil)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
         
        
        