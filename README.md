mersenne-twister-rs
==============

Pseudo random number generator [1] based on the [64-bit Mersenne Twister][2].
Rust port from C original [3][4].

References
----------
[1]: [The Art of Computer Programming, Donald E. Knuth, Vol. 2, Chapter 3](https://en.wikipedia.org/wiki/The_Art_of_Computer_Programming) <br/>
[2]: [Mersenne Twister] (https://en.wikipedia.org/wiki/Mersenne_Twister) <br/>
[3]: [Tables of 64-bit Mersenne Twisters by TAKUJI NISHIMURA](https://dl.acm.org/doi/pdf/10.1145/369534.369540) <br/>
[4]: [C code](http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/C-LANG/mt19937-64.c)


Run
-----

```
% cargo run 

 7266447313870364031
 4946485549665804864
16945909448695747420
16394063075524226720
 4873882236456199058
14877448043947020171
 6740343660852211943
13857871200353263164
 5249110015610582907
10205081126064480383
 1235879089597390050
17320312680810499042
16489141110565194782
 8942268601720066061
13520575722002588570
14226945236717732373
 9383926873555417063
15690281668532552105
11510704754157191257
15864264574919463609
 6489677788245343319
 5112602299894754389
10828930062652518694
15942305434158995996
15445717675088218264
 4764500002345775851
14673753115101942098
  236502320419669032
13670483975188204088
14931360615268175698
 8904234204977263924
12836915408046564963
12120302420213647524
15755110976537356441
 5405758943702519480
10951858968426898805
17251681303478610375
 4144140664012008120
18286145806977825275
13075804672185204371
10831805955733617705
 6172975950399619139
12837097014497293886
12903857913610213846
  560691676108914154
 1074659097419704618
14266121283820281686
11696403736022963346
13383246710985227247
 7132746073714321322
10608108217231874211
 9027884570906061560
12893913769120703138
15675160838921962454
 2511068401785704737
14483183001716371453
 3774730664208216065
 5083371700846102796
 9583498264570933637
17119870085051257224
 5217910858257235075
10612176809475689857
 1924700483125896976
 7171619684536160599
10949279256701751503
15596196964072664893
14097948002655599357
  615821766635933047
 5636498760852923045
17618792803942051220
  580805356741162327
  425267967796817241
 8381470634608387938
13212228678420887626
16993060308636741960
  957923366004347591
 6210242862396777185
 1012818702180800310
15299383925974515757
17501832009465945633
17453794942891241229
15807805462076484491
 8407189590930420827
  974125122787311712
 1861591264068118966
  997568339582634050
18046771844467391493
17981867688435687790
 3809841506498447207
 9460108917638135678
16172980638639374310
  958022432077424298
 4393365126459778813
13408683141069553686
13900005529547645957
15773550354402817866
16475327524349230602
 6260298154874769264
12224576659776460914
 6405294864092763507
 7585484664713203306
 5187641382818981381
12435998400285353380
13554353441017344755
  646091557254529188
11393747116974949255
16797249248413342857
15713519023537495495
12823504709579858843
 4738086532119935073
 4429068783387643752
  585582692562183870
 1048280754023674130
 6788940719869959076
11670856244972073775
 2488756775360218862
 2061695363573180185
 6884655301895085032
 3566345954323888697
12784319933059041817
 4772468691551857254
 6864898938209826895
 7198730565322227090
 2452224231472687253
13424792606032445807
10827695224855383989
11016608897122070904
14683280565151378358
 7077866519618824360
17487079941198422333
 3956319990205097495
 5804870313319323478
 8017203611194497730
 3310931575584983808
 5009341981771541845
 6930001938490791874
14415278059151389495
11001114762641844083
 6715939435439735925
  411419160297131328
 4522402260441335284
 3381955501804126859
15935778656111987797
 4345051260540166684
13978444093099579683
 9219789505504949817
 9245142924137529075
11628184459157386459
 7242398879359936370
 8511401943157540109
11948130810477009827
 6865450671488705049
13965005347172621081
15956599226522058336
 7737868921014130584
 2107342503741411693
15818996300425101108
16399939197527488760
13971145494081508107
 3910681448359868691
 4249175367970221090
 9735751321242454020
12418107929362160460
  241792245481991138
 5806488997649497146
10724207982663648949
 1121862814449214435
 1326996977123564236
 4902706567834759475
12782714623891689967
 7306216312942796257
15681656478863766664
  957364844878149318
 5651946387216554503
 8197027112357634782
 6302075516351125977
13454588464089597862
15638309200463515550
10116604639722073476
12052913535387714920
 2889379661594013754
15383926144832314187
 7841953313015471731
17310575136995821873
 9820021961316981626
15319619724109527290
15349724127275899898
10511508162402504492
 6289553862380300393
15046218882019267110
11772020174577005930
 3537640779967351792
 6801855569284252424
17687268231192623388
12968358613633237218
 1429775571144180123
10427377732172208413
12155566091986788996
16465954421598296115
12710429690464359999
 9547226351541565595
12156624891403410342
 2985938688676214686
18066917785985010959
 5975570403614438776
11541343163022500560
11115388652389704592
 9499328389494710074
 9247163036769651820
 3688303938005101774
 2210483654336887556
15458161910089693228
 6558785204455557683
 1288373156735958118
18433986059948829624
 3435082195390932486
16822351800343061990
 3120532877336962310
16681785111062885568
 7835551710041302304
 2612798015018627203
15083279177152657491
 6591467229462292195
10592706450534565444
 7438147750787157163
  323186165595851698
 7444710627467609883
 8473714411329896576
 2782675857700189492
 3383567662400128329
 3200233909833521327
12897601280285604448
 3612068790453735040
 8324209243736219497
15789570356497723463
 1083312926512215996
 4797349136059339390
 5556729349871544986
18266943104929747076
 1620389818516182276
  172225355691600141
 3034352936522087096
 1266779576738385285
 3906668377244742888
 6961783143042492788
17159706887321247572
 4676208075243319061
10315634697142985816
13435140047933251189
  716076639492622016
13847954035438697558
 7195811275139178570
10815312636510328870
 6214164734784158515
16412194511839921544
 3862249798930641332
 1005482699535576005
 4644542796609371301
17600091057367987283
 4209958422564632034
 5419285945389823940
11453701547564354601
 9951588026679380114
 7425168333159839689
 8436306210125134906
11216615872596820107
 3681345096403933680
 5770016989916553752
11102855936150871733
11187980892339693935
  396336430216428875
 6384853777489155236
 7551613839184151117
16527062023276943109
13429850429024956898
 9901753960477271766
 9731501992702612259
 5217575797614661659
10311708346636548706
15111747519735330483
 4353415295139137513
 1845293119018433391
11952006873430493561
 3531972641585683893
16852246477648409827
15956854822143321380
12314609993579474774
16763911684844598963
16392145690385382634
 1545507136970403756
17771199061862790062
12121348462972638971
12613068545148305776
  954203144844315208
 1257976447679270605
 3664184785462160180
 2747964788443845091
15895917007470512307
15552935765724302120
16366915862261682626
 8385468783684865323
10745343827145102946
 2485742734157099909
  916246281077683950
15214206653637466707
12895483149474345798
 1079510114301747843
10718876134480663664
 1259990987526807294
 8326303777037206221
14104661172014248293
15531278677382192198
 3874303698666230242
 3611366553819264523
 1358753803061653874
 1552102816982246938
14492630642488100979
15001394966632908727
 2273140352787320862
17843678642369606172
 2903980458593894032
16971437123015263604
12969653681729206264
 3593636458822318001
 9719758956915223015
 7437601263394568346
 3327758049015164431
17851524109089292731
14769614194455139039
 8017093497335662337
12026985381690317404
  739616144640253634
15535375191850690266
 2418267053891303448
15314073759564095878
10333316143274529509
16565481511572123421
16317667579273275294
13991958187675987741
 3753596784796798785
 9078249094693663275
 8459506356724650587
12579909555010529099
 7827737296967050903
 5489801927693999341
10995988997350541459
14721747867313883304
 7915884580303296560
 4105766302083365910
12455549072515054554
13602111324515032467
 5205971628932290989
 5034622965420036444
 9134927878875794005
11319873529597990213
14815445109496752058
 2266601052460299470
 5696993487088103383
 6540200741841280242
 6631495948031875490
 5328340585170897740
17897267040961463930
 9030000260502624168
14285709137129830926
12854071997824681544
15408328651008978682
 1063314403033437073
13765209628446252802
  242013711116865605
 4772374239432528212
 2515855479965038648
 5872624715703151235
14237704570091006662
  678604024776645862
12329607334079533339
17570877682732917020
 2695443415284373666
 4312672841405514468
 6454343485137106900
 8425658828390111343
16335501385875554899
 5551095603809016713
11781094401885925035
 9395557946368382509
 9765123360948816956
18107191819981188154
16049267500594757404
16349966108299794199
 1040405303135858246
 2366386386131378192
  223761048139910454
15375217587047847934
15231693398695187454
12916726640254571028
 8878036960829635584
 1626201782473074365
 5758998126998248293
18077917959300292758
10585588923088536745
15072345664541731497
 3559348759319842667
12744591691872202375
 2388494115860283059
 6414691845696331748
 3069528498807764495
 8737958486926519702
18059264986425101074
 3139684427605102737
12378931902986734693
  410666675039477949
12139894855769838924
 5780722552400398675
 7039346665375142557
 3020733445712569008
 2612305843503943561
13651771214166527665
16478681918975800939
  566088527565499576
 4715785502295754870
 6957318344287196220
11645756868405128885
13139951104358618000
17650948583490040612
18168787973649736637
 5486282999836125542
 6122201977153895166
17324241605502052782
10063523107521105867
17537430712468011382
10828407533637104262
10294139354198325113
12557151830240236401
16673044307512640231
10918020421896090419
11077531235278014145
 5499571814940871256
 2334252435740638702
18177461912527387031
 2000007376901262542
 7968425560071444214
 1472650787501520648
 3115849849651526279
 7980970700139577536
12153253535907642097
 8109716914843248719
 3154976533165008908
 5553369513523832559
10345792701798576501
 3677445364544507875
10637177623943913351
 7380255087060498096
14479400372337014801
15381362583330700960
  204531043189704802
13699106540959723942
 3817903465872254783
10972364467110284934
 2701394334530963810
 2931625600749229147
16428252083632828910
11873166501966812913
 5566810080537233762
 7840617383807795056
10699413880206684652
18259119259617231436
10332714341486317526
10137911902863059694
  669146221352346842
 8373571610024623455
10620002450820868661
12220730820779815970
 5902974968095412898
 7931010481705150841
16413777368097063650
11273457888324769727
13719113891065284171
 8327795098009702553
10333342364827584837
 6202832891413866653
 9137034567886143162
14514450826524340059
  473610156015331016
  813689571029117640
13776316799690285717
10429708855338427756
 8995290140880620858
 2320123852041754384
 8082864073645003641
 6961777411740398590
10008644283003991179
 3239064015890722333
16762634970725218787
16467281536733948427
10563290046315192938
 5108560603794851559
15121667220761532906
14155440077372845941
10050536352394623377
15474881667376037792
 3448088038819200619
 3692020001240358871
 6444847992258394902
 8687650838094264665
 3028124591188972359
16945232313401161629
15547830510283682816
 3982930188609442149
14270781928849894661
13768475593433447867
13815150225221307677
 8502397232429564693
  718377350715476994
 7459266877697905475
 8353375565171101521
 7807281661994435472
16924127046922196149
10157812396471387805
 2519858716882670232
 7384148884750265792
 8077153156180046901
 3499231286164597752
 2700106282881469611
14679824700835879737
14188324938219126828
 3016120398601032793
10858152824243889420
 9412371965669250534
 4857522662584941069
  984331743838900386
 4094160040294753142
 2368635764350388458
15101240511397838657
15584415763303953578
 7831857200208015446
 1952643641639729063
 4184323302594028609
16795120381104846695
 3541559381538365280
15408472870896842474
 5628362450757896366
16277348886873708846
12437047172652330846
10172715019035948149
 1999700669649752791
 6217957085626135027
11220551167830336823
16478747645632411810
 5437280487207382147
11382378739613087836
15866932785489521505
 5502694314775516684
16440179278067648435
15510104554374162846
15722061259110909195
10760687291786964354
10736868329920212671
 4166148127664495614
14303518358120527892
 9122250801678898571
10028508179936801946
  216630713752669403
10655207865433859491
 4041437116174699233
 6280982262534375348
  297501356638818866
13976146806363377485
13752396481560145603
11472199956603637419
16393728429143900496
14752844047515986640
 1524477318846038424
 6596889774254235440
 1591982099532234960
 8065146456116391065
 3964696017750868345
17040425970526664920
11511165586176539991
 3443401252003315103
16314977947073778249
16860120454903458341
 5370503221561340846
15362920279125264094
 2822458124714999779
14575378304387898337
 9689406052675046032
 2872149351415175149
13019620945255883050
14929026760148695825
 8503417349692327218
 9677798905341573754
  828949921821462483
16110482368362750196
15794218816553655671
14942910774764855088
12026350906243760195
13610867176871462505
18324536557697872582
 2658962269666727629
  327225403251576027
 9207535177029277544
 8744129291351887858
 6129603385168921503
18385497655031085907
13024478718952333892
14547683159720717167
 5932119629366981711
  325385464632594563
 3559879386019806291
 6629264948665231298
14358245326238118181
15662449672706340765
13975503159145803297
 3609534220891499022
 4224273587485638227
 9274084767162416370
13156843921244091998
18284750575626858789
14664767920489118779
11292057742031803221
13919998707305829132
14473305049457001422
 9696877879685767807
 1406758246007973837
 2429517644459056881
14361215588101587430
11386164476149757528
10474116023593331839
 2921165656527786564
15604610369733358953
12955027028676000544
10314281035410779907
 3167047178514709947
 1088721329408346700
17930425515478182741
 7466411836095405617
15534027454610690575
10879629128927506091
11502219301371200635
13915106894453889418
 4226784327815861027
12335222183627106346
 3648499746356007767
18441388887898023393
18117929843327093625
 4237736098094830438
14229123019768296655
 3930112058127932690
12663879236019645778
 9281161952002617309
 4978473890680876319
  845759387067546611
 1386164484606776333
 8008554770639925512
11159581016793288971
18065390393740782906
17647985458967631018
 9092379465737744314
 2914678236848656327
 4376066698447630270
16057186499919087528
 3031333261848790078
 2926746602873431597
 7931945763526885287
  147649915388326849
15801792398814946230
 5265900391686545347
16173686275871890830
 7562781050481886043
 5853506575839330404
14957980734704564792
10944286556353523404
 1783009880614150597
 9529762028588888983
  822992871011696119
 2130074274744257510
 8000279549284809219
 3514744284158856431
  128770032569293263
 3737367602618100572
16364836605077998543
  783266423471782696
 4569418252658970391
11093950688157406886
14888808512267628166
 4217786261273670948
17047486076688645713
14133826721458860485
17539744882220127106
12394675039129853905
 5757634999463277090
 9621947619435861331
 1182210208559436772
14603391040490913939
17481976703660945893
14063388816234683976
 2046622692581829572
 8294969799792017441
 5293778434844788058
17976364049306763808
  399482430848083948
16495545010129798933
15241340958282367519
  989828753826900814
17616558773874893537
 2471817920909589004
11764082277667899978
 9618755269550400950
 1240014743757147125
 1887649378641563002
 1842982574728131416
13243531042427194002
 7688268125537013927
 3080422097287486736
 2562894809975407783
12428984115620094788
 1355581933694478148
 9895969242586224966
 8628445623963160889
 4298916726468199239
12773165416305557280
 5240726258301567487
 4975412836403427561
 1842172398579595303
 7812151462958058676
17974510987263071769
14980707022065991200
18294903201142729875
12911672684850242753
 8979482998667235743
16808468362384462073
 5981317232108359798
12373702800369335100
16119707581920094765
 2782738549717633602
15454155188515389391
16495638000603654629
16348757069342790497
 7769562861984504567
17504300515449231559
 5557710032938318996
11846125204788401203
13957316349928882624
 2738350683717432043
15738068448047700954
 6224714837294524999
 6081930777706411111
11366312928059597928
 4355315799925031482
12393324728734964015
15277140291994338591
 1406052433297386355
15859448364509213398
 1672805458341158435
 2926095111610982994
11056431822276774455
12083767323511977430
 3296968762229741153
12312076899982286460
17769284994682227273
15349428916826953443
 1056147296359223910
18305757538706977431
 6214378374180465222
14279648441175008454
17791306410319136644
  956593013486324072
 2921235772936241950
10002890515925652606
10399654693663712506
 6446247931049971441
 6380465770144534958
11439178472613251620
10131486500045494660
 3692642123868351947
10972816599561388940
 4931112976348785580
 8213967169213816566
15336469859637867841
15026830342847689383
 7524668622380765825
17309937346758783807
  372780684412666438
 5642417144539399955
18303842993081194577
11085303253831702827
15658163165983586950
 8517521928922081563
16091186344159989860
17614656488010863910
 4736067146481515156
13449945221374241354
17755469346196579408
13300502638545717375
 6611828134763118043
14177591906740276597
 9340430243077460347
 7499765399826404087
 3409518087967832469
 9013253864026602045
 4444307427984430192
 3729283608700519712
13642048880719588383
16486557958022946240
 2996465014991157904
10020049344596426576
12302485648009883778
 8492591321344423126
17407986443716172520
10530482934957373052
15740662350540828750
 1790629986901049436
 6305948377669917188
15092985352503125323
  928505047232899787
14404651977039851607
 7564177565277805597
 3411236815351677870
 7752718145953236134
12315979971311483798
12477729506691004724
14654956300924793305
 6689803038918974388
 1540738812233000153
13508351811701989957
15864432023192136053
 7990997967273843917
 7424300239290765161
   39585249496300263
 3877436595063283319
10710642254398044448
 4653804418844456375
 1232267496410380283
 3690525514009038824
15459770765077428485
13240346522153894145
 5674964360688390624
16973644653010587289
15924280764204855206
15196708627253442662
17596174821341373274
16196745023027393691
 6980050627399795351
17582264380857746637
18170372407506856324
12108126025631005514
15687749089493373169
 5814107289258228434
 9381977959648494876
15895601183088112734
16267869075651604263
15228381979765852785
11949618678312581999
 4545324791131029438
  582725409406225185
15282520250746126790
14758446535973412711
 7605613563088071833
 1111140641057375915
 5364843095234852245
  218335432181198977
 4891472444796201742
 4564628942836375772
15500501278323817088
 4913946328556108657
 2684786251736694229
12090498456116310122
 5310885782157038567
 5032788439854011923
12627401038822728242
11869662610126430929
17650156853043540226
12126672500118808436
10437658933435653256
13133995470637873311
 4601324715591152820
 1874350460376708372
 5808688626286061164
13777088437302430376
 5018451954762213522
 2588296738534474754
 5503414509154170711
 5230497186769951796
13261090710400573914
 8515217303152165705
11074538219737365303
15481562385740613213
12705484409881007350
14221931471178549498
12905633420087112297
17337759164357146506
14081997515778175224
17384320185513122939
 7131793076779216692
17483217190312403109
  900692047897995877
14723287313048560400
 6132094372965340305
 7572797575350925726
12725160700431903514
  380860122911632449
 1900504978569024571
 8423729759529914138
 7305587201606052334
12446871355267313320
 4615812356515386206
 3361817115406652303
17690418922000878428
14632214537567910559
 2709702289926174775
 3459675155951086144
 7788364399926538150
16043992474431955950
15830963823784930267
 4216893617835797954
  538159724689093771
16029152738918251363
14444848757576686696
12941757045272633696
10900480525147953314
12547307449905859302
16001571796892398181
  407942194622690676
13873235372903944444
18071603799493008777
 1015646077646778622
 9387605808959554815
11566702442022019410
 7061722181092883183
 2629032108249254109
 5271820053177594520
12640880742139693547
10098688629735675775
 5716304472850923064
 3312674502353063071
 7295926377425759633
  833281439103466115
16316743519466861667
 9912050326606348167
11651133878100804242
18026798122431692459
 6157758321723692663
 4856021830695749349
 7074321707293278978
10748097797809573561
 2949954440753264783
 9813922580940661152
 9949237950172138336
15643982711269455885
16078663425810239127
12508044395364228880
12920301578340189344
15368071871011048915
 1610400750626363239
11994736084146033126
 6042574085746186088
 4154587549267685807
15915752367312946034
 1191196620621769193
  467437822242538360
 2836463788873877488
10476401302029164984
 1716169985450737419
 5327734953288310341
 3994170067185955262
  884431883768190063
11019001754831208284
14322807384384895215
  161011537360955545
 1466223959660131656
 5227048585229497539
12410731857504225031
 2142243279080761103
17682826799106851430
 1792612570704179953
14727410295243056025
 1459567192481221274
 5669760721687603135
17507918443756456845
10354471145847018200
10362475129248202288
13143844410150939443
 6861184673150072028
18396524361124732580
  543906666394301875
12476817828199026728
11853496871128122868
12747674713108891748
 7986179867749890282
 9158195177777627533
 2217320706811118570
 8631389005200569973
 5538133061362648855
 3369942850878700758
 7813559982698427184
  509051590411815948
10197035660403006684
13004818533162292132
 9831652587047067687
 7619315254749630976
  994412663058993407
0.35252031
0.51052342
0.79771733
0.39300273
0.27216673
0.72151068
0.43144703
0.38522290
0.20270676
0.58227313
0.80812143
0.83767297
0.92401619
0.84065425
0.00852052
0.13975395
0.35250930
0.71196972
0.14627395
0.17775331
0.61046382
0.49623272
0.23292425
0.25038837
0.04380664
0.43275994
0.74540936
0.33830700
0.68832616
0.68744230
0.63626548
0.85932936
0.37089670
0.50756304
0.69925960
0.83481025
0.09053196
0.09523253
0.17783108
0.78027239
0.70071054
0.51879252
0.83027285
0.92895011
0.72144803
0.18868644
0.83655674
0.20358945
0.99852143
0.88340103
0.46729949
0.96993433
0.00162682
0.46829774
0.59080423
0.54921999
0.42516462
0.54952196
0.99534722
0.04473888
0.71139235
0.91881407
0.33781561
0.45746234
0.78292126
0.69206723
0.66175448
0.07091147
0.18179208
0.38168454
0.38819527
0.42452711
0.22732724
0.16191307
0.36842667
0.13060083
0.68833248
0.60498705
0.19195304
0.26628584
0.17030858
0.23892426
0.38430236
0.28034283
0.76069020
0.21560653
0.78101667
0.90847812
0.06467974
0.18487868
0.23570471
0.29475460
0.65563767
0.10066446
0.57272419
0.88731391
0.60650995
0.96346079
0.32940100
0.29977746
0.03798193
0.18026822
0.22402746
0.45480119
0.98114604
0.25800668
0.94362433
0.17901062
0.36019313
0.45933644
0.68309457
0.28175454
0.00774729
0.77054527
0.99723413
0.59807532
0.10294164
0.32429228
0.54928986
0.18410980
0.08441555
0.14230333
0.58892064
0.94030475
0.35378784
0.77584320
0.71222448
0.83565208
0.47309248
0.23810761
0.74408520
0.08891527
0.09729786
0.38377368
0.05092308
0.69065638
0.10449489
0.45050670
0.92209534
0.80083714
0.27902692
0.26897142
0.50650468
0.80111472
0.54590012
0.96406097
0.63779553
0.81054357
0.75369248
0.47473037
0.89100315
0.89395984
0.09985519
0.34087631
0.22293557
0.24375510
0.31764191
0.04076993
0.06160830
0.41333434
0.11883030
0.04548820
0.01008040
0.25336184
0.07325432
0.49860151
0.07148695
0.89483338
0.87054457
0.15116809
0.59650469
0.47487776
0.43490298
0.36684681
0.16470796
0.76865078
0.42920071
0.20545481
0.87615922
0.80332404
0.36462506
0.49571309
0.51904488
0.15534589
0.43719893
0.16562157
0.37290862
0.91842631
0.21310523
0.87849154
0.18532269
0.81713354
0.52182344
0.51845619
0.96261204
0.18758718
0.68897600
0.61484764
0.46752993
0.05865458
0.11614359
0.90386866
0.45781805
0.70649579
0.50917048
0.21210656
0.97818608
0.00788342
0.61375222
0.67366318
0.24197878
0.66177985
0.10463932
0.67390799
0.50025262
0.88332650
0.77966851
0.13403622
0.54357114
0.97664854
0.06540961
0.24013176
0.67234032
0.91347883
0.35486839
0.87207865
0.43036581
0.23652488
0.81238450
0.72058432
0.42239916
0.80265764
0.03552838
0.61939480
0.50972420
0.21053832
0.59952743
0.36821802
0.45659617
0.12529468
0.76941623
0.99878168
0.08602783
0.81825937
0.39350710
0.86090923
0.36090230
0.75628888
0.45036982
0.44602266
0.20595631
0.62241953
0.36777732
0.47523727
0.50248178
0.73570362
0.48237781
0.45590948
0.73580783
0.96403851
0.94586342
0.48819868
0.48102038
0.94618182
0.90279924
0.78396650
0.85182389
0.92149394
0.32679198
0.83554856
0.28320609
0.34598409
0.82090005
0.40177958
0.38888785
0.77873931
0.23297931
0.75329335
0.30770340
0.71417540
0.68939065
0.36577776
0.50784857
0.50928090
0.02552055
0.85999075
0.26692089
0.01402799
0.67550392
0.48305605
0.74608351
0.63408891
0.58904230
0.44337996
0.42174728
0.74041679
0.72719148
0.19801992
0.66263633
0.10381594
0.32818760
0.68369661
0.56076212
0.68681921
0.91616269
0.39836106
0.39685027
0.97507945
0.91010563
0.27447360
0.95538357
0.76758522
0.60091060
0.37734461
0.82948248
0.06598078
0.50147615
0.08417763
0.18910044
0.51661735
0.55011011
0.64888175
0.82986845
0.15126656
0.92649390
0.25494941
0.73275293
0.94184393
0.84755226
0.45921936
0.72934054
0.43722403
0.34305596
0.10827860
0.29026676
0.01935431
0.46668573
0.83247509
0.26349603
0.01938542
0.43222250
0.18109983
0.29337450
0.16721917
0.94751650
0.67795254
0.56666228
0.20699452
0.23247262
0.19138610
0.73495506
0.85893600
0.83411526
0.93689655
0.91804752
0.99352333
0.03207550
0.28386071
0.48029543
0.18736013
0.31736452
0.72542230
0.57530912
0.04229918
0.84798296
0.21886935
0.98655615
0.52243102
0.22611020
0.42975741
0.21726739
0.10912048
0.96684473
0.01092456
0.12461901
0.57989070
0.39848707
0.06330277
0.62826828
0.01159081
0.23157320
0.64690912
0.44876902
0.04463930
0.18933780
0.21284518
0.61363480
0.67144845
0.38625586
0.75719122
0.40361050
0.26708873
0.54534727
0.90174015
0.58654140
0.44885346
0.35505544
0.65317830
0.26074572
0.39472912
0.54366914
0.75020660
0.76113614
0.24595582
0.03941247
0.60356153
0.23615721
0.01603475
0.72432457
0.39837424
0.04195329
0.81561058
0.34208440
0.00513953
0.92826234
0.11410393
0.86692030
0.25238726
0.98258626
0.53353856
0.72269001
0.71850984
0.66829681
0.03540769
0.01676450
0.23557835
0.78758497
0.85969589
0.14673207
0.28013860
0.17796942
0.69924087
0.44663597
0.62112513
0.44079883
0.48995231
0.18411497
0.18440877
0.74016388
0.28845694
0.22969080
0.76851164
0.15551473
0.28980810
0.40906710
0.47619039
0.72611392
0.55802939
0.69365597
0.85736313
0.83343150
0.21324760
0.45327806
0.33053855
0.98198279
0.53279389
0.76877035
0.20548656
0.37065042
0.59026910
0.67418036
0.23585843
0.98156397
0.27849804
0.56198954
0.68752287
0.30073445
0.69348664
0.72515585
0.40629047
0.09320027
0.24334978
0.91407662
0.97226538
0.33904970
0.01717092
0.60155725
0.03001652
0.50979706
0.80531036
0.17450719
0.84984399
0.00498130
0.51636405
0.14080868
0.62289701
0.07853030
0.70567541
0.79844050
0.63766566
0.03559031
0.40994535
0.08423996
0.00389626
0.50608347
0.19622681
0.90537903
0.75458034
0.75102094
0.81491673
0.92925931
0.38074332
0.54817053
0.72593246
0.02146791
0.57990460
0.87921074
0.59913886
0.66726893
0.24269154
0.73344575
0.71826052
0.92313935
0.05212996
0.93771536
0.69489385
0.57581887
0.48106155
0.06808800
0.33633940
0.69142320
0.46566781
0.70654143
0.16541368
0.76257631
0.82777900
0.62958327
0.34757935
0.10891487
0.79912728
0.01156543
0.23111261
0.58535640
0.87461956
0.21723454
0.80409615
0.33169686
0.72800785
0.31218099
0.13729737
0.41637635
0.01234597
0.58313811
0.66746028
0.05105595
0.14930937
0.56044864
0.76196851
0.98800104
0.37075949
0.88740864
0.40697115
0.96598278
0.86013661
0.85386784
0.23986516
0.39027464
0.59593927
0.00161530
0.31768197
0.65702729
0.66461914
0.62937471
0.92120758
0.87578958
0.37539860
0.59182615
0.12092214
0.55130437
0.86365117
0.38725162
0.28757657
0.42803199
0.39014405
0.50253853
0.85306128
0.92018995
0.71421618
0.54236780
0.96221157
0.22956898
0.96519876
0.06694102
0.11915854
0.01354308
0.24720070
0.71671739
0.00604305
0.65012352
0.71151390
0.46616159
0.99228224
0.20684576
0.62941006
0.84535326
0.30678993
0.55264568
0.50094784
0.39409122
0.15479416
0.36536318
0.51925656
0.65567178
0.67255519
0.55089659
0.42194295
0.27172413
0.79540954
0.71594806
0.88372598
0.29179452
0.66411306
0.57064687
0.42494633
0.73389255
0.12097313
0.53338622
0.38493233
0.79348021
0.01851341
0.58594454
0.88396240
0.04410730
0.67419924
0.62770011
0.64644200
0.40335135
0.17952644
0.55564678
0.56643922
0.37715015
0.87092180
0.56726159
0.34011210
0.13661819
0.11474177
0.93930097
0.48549077
0.28484289
0.13374371
0.40966056
0.73662873
0.37355323
0.65216092
0.27372469
0.56032082
0.14882684
0.95462890
0.17090266
0.92374766
0.98368259
0.68448367
0.02872548
0.68598279
0.04601084
0.17170501
0.08906644
0.23730372
0.02929037
0.38566261
0.68957569
0.53021050
0.44200157
0.32085701
0.72520053
0.17454174
0.19676599
0.88243877
0.87030228
0.15124486
0.78670160
0.51731632
0.56674531
0.20910664
0.84962640
0.05220467
0.91783159
0.19138968
0.68126378
0.79574471
0.14910848
0.28030331
0.98067264
0.31263980
0.67448964
0.69266650
0.40033551
0.22789781
0.78317066
0.55815261
0.11247054
0.47337901
0.46310033
0.53192452
0.56164078
0.41750378
0.43880622
0.69739327
0.11092778
0.18333765
0.67222441
0.12789170
0.88316806
0.37891271
0.14935268
0.64522185
0.93902079
0.62481092
0.21794927
0.71535266
0.62169579
0.65147153
0.01411645
0.96413465
0.01021578
0.50605180
0.51595053
0.03308040
0.01497870
0.07809658
0.35743383
0.58079701
0.11785557
0.89568677
0.38793964
0.37117709
0.13994133
0.11032813
0.99998594
0.06695042
0.79774786
0.11093584
0.23879095
0.85918615
0.16109636
0.63479696
0.75023359
0.29061187
0.53764772
0.30652318
0.51387302
0.81620973
0.82433610
0.18302488
0.79048957
0.07598187
0.27887732
0.37061042
0.36441016
0.93736882
0.77480946
0.02269132
0.40309874
0.16427650
0.13969296
0.57605029
0.00242426
0.56626691
0.84390990
0.87455806
0.12321023
0.87561663
0.60431578
0.35880839
0.50426282
0.50697689
0.06631164
0.14976092
0.89356018
0.91473662
0.04235237
0.50073724
0.75969690
0.91743994
0.79352335
0.58078351
0.91819984
0.53520520
0.18267367
0.05608828
0.68315721
0.27264599
0.41245634
0.69706222
0.69666203
0.08967342
0.64081905
0.22576796
0.69315628
0.53981640
0.76059129
0.56712344
0.94318621
0.44081094
0.31699284
0.29477911
0.80069824
0.28366921
0.96718081
0.85345644
0.11681215
0.47600710
0.33448255
0.31217271
0.35469241
0.59511650
0.49583692
0.48922303
0.20215259
0.60159380
0.17882055
0.77601258
0.71020391
0.41833503
0.71522856
0.87534517
0.43703394
0.43056077
0.64828071
0.43069441
0.39356849
0.32063367
0.92788963
0.16878266
0.56762591
0.56042446
0.84958464
0.79408949
0.08220340
0.13922856
0.82529019
0.27134959
0.00278080
0.66192389
0.01782933
0.95404763
0.50787645
0.85320521
0.83690362
0.83771227
0.46268665
0.31716742
0.01716647
0.68264674
0.01789888
0.30446846
0.14942271
0.26982182
0.74933947
0.50394161
0.78444542
0.40009256
0.40333422
0.16627342
0.01898760
0.04221829
0.77960213
0.66230976
0.56015996
0.49535426
0.38536259
0.40406773
0.99930568
0.00857945
0.16158390
0.64805163
0.20237524
0.59106326
0.76968277
0.96887042
0.29264851
0.97373775
0.16767633
0.33014482
0.27426548
0.10947014
0.75920652
0.37757457
0.13125207
0.00826451
0.96684342
0.69362226
0.22763554
0.20717541
0.42112268
0.22803038
0.33481806
0.14968742
0.71598558
0.55126711
0.64518015
0.65170197
0.89103003
0.72728361
0.24485454
0.09410780
0.79818029
0.54212409
0.17790462
0.64442619
0.62193511
0.51193256
0.02848781
0.05719604
0.45795152
0.03219332
0.28310254
0.85746127
0.64890240
0.20658356
0.50946422
0.80432490
0.08354468
0.09222723
0.67455943
0.44638771
0.76366629
0.99677267
0.89311242
0.11627279
0.09181302
0.44767077
0.16448724
0.26005539
0.28670391
0.52465703
0.43598116
0.41869096
0.98043420
0.01497272
0.51791571
0.61825308
0.85503436
0.63025655
0.02719292
0.09865668
0.30321729
0.56998039
0.14946350
0.64823918
0.19931639
0.14623555
0.54169913
0.68944135
0.73551005
0.46743658
0.04109096
0.26625801
0.09537298
0.98207890
0.58109721
0.70793680
0.84379365
0.42774726
0.12653597
0.08566633
0.53366781
0.33960092
0.11036831
0.84464510
0.16493476
0.92493443
0.87640673
0.52727644
0.57181349
0.65071340
0.00978637
0.31700693
0.69148222
0.85063311
0.06781819
0.30794534
0.65541667
0.16400484
0.06886223
0.96227205
0.09633060
0.34513153
0.31013900
0.78165882
0.39583699
0.86327936
0.69269199
0.11016575
0.67358419
0.81775427
0.50052824
0.30068582
0.16606837
0.62243724
0.47863741
0.68796498
0.31526949
0.41180883
0.23022147
0.82342139
0.83003381
0.53571829
0.41081533
0.48600142
```

