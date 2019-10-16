load(file = "testdata.rda")

### Test runSingleTraitGwas

stg0 <- runSingleTraitGwas(gData = gDataTest, trials = 1)
stg1 <- runSingleTraitGwas(gData = gDataTest)
result1 <- runSingleTraitGwas(gData = gDataTest, trials = 1,
                              covar = "V1")[["GWAResult"]]
result2 <- runSingleTraitGwas(gData = gDataTest, trials = 1,
                              snpCov = "M2")[["GWAResult"]]
result3 <- runSingleTraitGwas(gData = gDataTest, trials = 1, covar = "V1",
                              snpCov = "M2")[["GWAResult"]]

## Check output structure.

expect_true(inherits(stg0, "GWAS"))
expect_equal(length(stg0), 5)
expect_equal(names(stg0), 
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal(names(stg0[["GWASInfo"]]), 
             c("call", "remlAlgo", "thrType", "MAF", "GLSMethod", "varComp", 
               "genomicControl", "inflationFactor"))
expect_equal(length(stg0[["GWASInfo"]][["varComp"]][["ph1"]]), 5)
expect_true(inherits(stg0[["GWAResult"]], "list"))
expect_equal(length(stg0[["GWAResult"]]), 1)
expect_equal(names(stg0[["GWAResult"]]), "ph1")
expect_equal(length(stg1[["GWAResult"]]), 2)
expect_equal(names(stg1[["GWAResult"]]), c("ph1", "ph2"))

## Check p-Values

expect_equal(stg0[["GWAResult"]][["ph1"]][["pValue"]],
             c(0.556410736473299, 0.644300498575488, 0.815015376305757, 
               0.46416067677484, 0.307927840927776, 0.937588447652476,
               0.854285085985704, 0.640313489847305, 0.489678729290149, 
               0.673232138421906, 0.552144063962775, 0.395714586737013, 
               0.832011110721433, 0.156889391235401, 0.253195176015333))
expect_equal(result1[["ph1"]][["pValue"]],
             c(0.351340586556164, 0.805121913672965, 0.97710572208884, 
               0.465199071828829, 0.32068791201796, 0.93835429411708,
               0.705896462282515, 0.728326844507144, 0.594118369168397,
               0.771670607342803, 0.489323081632407, 0.473509357066004, 
               0.871195025941687, 0.151626281146346, 0.233293824250749))
expect_equal(result2[["ph1"]][["pValue"]],
             c(0.601870706247044, 0.639608839221507, 0.879702684738528,
               0.528671908918732, 0.307927840927776, 0.805725551506335,
               0.827794265956364, 0.620674116521037, 0.579960510641916, 
               0.729405950213325, 0.552144063962774, 0.354221370612768, 
               0.934748044050838, 0.156889391235401, 0.328324199598967))
expect_equal(result3[["ph1"]][["pValue"]],
             c(0.407453651106354, 0.805121913672964, 0.964222347313319, 
               0.605576314394518, 0.320687912017961, 0.878249328147975, 
               0.735648828647385, 0.715294618992726, 0.652238864565491,
               0.893223912521159, 0.489323081632408, 0.457210861098801, 
               0.880519035238243, 0.151626281146346, 0.239301792260896))

## Check effects.

expect_equal(stg0[["GWAResult"]][["ph1"]][["effect"]],
             c(0.281293915064906, 0.270883214623018, -0.183866889219171, 
               -0.424623375012714, -0.712404445545492, -0.0754950435841936, 
               -0.127726819554899, -0.389474297834824, 0.758085652353884,
               0.177147797855602, 0.30408544635742, 0.572138780492293,
               -0.156479414190887, -1.19930737632641, 1.32017648576984))
expect_equal(result1[["ph1"]][["effect"]],
             c(0.478643894550101, 0.156154597676827, 0.0247115254274738, 
               -0.484776962765978, -0.774132902165756, -0.0849665260154638, 
               -0.297031691751697, -0.321370575524584, 0.664037543019707,
               0.138940358049405, 0.389553570358166, 0.547438870570076,
               -0.137123618369648, -1.34196336062861, 1.55024604899299))
expect_equal(result2[["ph1"]][["effect"]],
             c(0.265287119537725, 0.274992779365586, -0.128016335697819,
               -0.369104717582792, -0.71240444554548, -0.240481834434491,
               -0.159371598035534, -0.415486560755866, 0.661357311094297,
               0.1534068496046, 0.304085446357424, 0.657987591567463, 
               -0.0571392499051401, -1.19930737632643, 1.08377527208679))
expect_equal(result3[["ph1"]][["effect"]],
             c(0.473146507899702, 0.15615459767689, 0.0419019448563874,
               -0.357345442796037, -0.774132902165716, -0.169946796542321,
               -0.290281388459914, -0.339538772577021, 0.611810266543856,
               0.0692437953496002, 0.38955357035814, 0.593439913315101, 
               0.121092108320429, -1.34196336062862, 1.4138231687998))

## Check values for traits containing NAs.

expect_equal(stg1[["GWAResult"]][["ph2"]][["pValue"]],
             c(0.535866686599578, 0.738175077654468, NA, 0.862874874084195, 
               0.707996660593566, 0.814639480900727, 0.938915044727933,
               0.704903948547425, 0.447048189560569, 0.894606032126416,
               0.511413255485449, 0.457037891963561, 0.966898351313243,
               0.253600829295148, 0.782993491205565))
expect_equal(stg1[["GWAResult"]][["ph2"]][["effect"]],
             c(0.384318156764552, 0.234891575473521, NA, 0.0992876321027796, 
               0.313594562157908, -0.197059831391177, -0.0582949760905214, 
               -0.329230600354315, 0.846952100377062, 0.0555585062439554,
               -0.3325198897839, 0.285348088200379, -0.0221226243845664,
               -0.810848162781112, -0.119937375119302))