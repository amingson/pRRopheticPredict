#在cpg2016 tissuetype中添加了allsolidtumor 这一选项，考虑当时作者忘记写上了
getCGPinfo.my <- function (drug, tissueType = "all", dataset = "cgp2014") 
{
  possibleDatasets <- c("cgp2014", "cgp2016")
  if (!dataset %in% possibleDatasets) 
    stop(paste("ERROR: the dataset specified was not found. Note dataset names are case sensitive. Please select from: ", 
               (paste(possibleDatasets, collapse = ", "))))
  if (dataset == "cgp2014") {
    if (!tissueType %in% c("all", "allSolidTumors", "blood", 
                           "breast", "CNS", "GI tract", "lung", "skin", "upper aerodigestive")) 
      stop("ERROR: the tissue type specified must be one of \"all\", \"allSolidTumors\", \"blood\", \"breast\", \"CNS\", \"GI tract\", \"lung\", \"skin\" or \"upper aerodigestive\". These tissue types are represented by at least 40 cell lines (although not all 40 may have been screened with each drug).")
    possibleDrugs <- c("A.443654", "A.770041", "ABT.263", 
                       "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", 
                       "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", 
                       "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", 
                       "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", 
                       "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", 
                       "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", 
                       "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", 
                       "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", 
                       "CCT007093", "CCT018159", "CEP.701", "CGP.082996", 
                       "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", 
                       "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", 
                       "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", 
                       "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", 
                       "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", 
                       "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", 
                       "GSK.650394", "GW.441756", "GW843682X", "Imatinib", 
                       "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", 
                       "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", 
                       "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", 
                       "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", 
                       "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", 
                       "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", 
                       "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", 
                       "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", 
                       "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", 
                       "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", 
                       "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", 
                       "SB.216763", "SB590885", "Shikonin", "SL.0101.1", 
                       "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", 
                       "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", 
                       "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", 
                       "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", 
                       "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", 
                       "ZM.447439")
    if (!drug %in% possibleDrugs) 
      stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", 
                 (paste(possibleDrugs, collapse = ", "))))
    data(drugAndPhenoCgp)
    colIc50Name <- paste(drug, "_IC_50", sep = "")
    ic50s <- as.numeric(drugSensitivityDataCgp[, colIc50Name])
    names(ic50s) <- drugSensitivityDataCgp[, "Cell.Line"]
    whichNas <- which(is.na(ic50s))
    ic50s <- ic50s[-whichNas]
    tissue <- drugSensitivityDataCgp[, "Cancer.Type"]
    names(tissue) <- drugSensitivityDataCgp[, "Cell.Line"]
    tissue <- tissue[-whichNas]
    if (tissueType != "all") {
      if (tissueType == "allSolidTumors") {
        tissueType <- ic50s <- ic50s[!(tissue %in% "blood")]
      }
      else {
        ic50s <- ic50s[tissue %in% tissueType]
      }
    }
    pDataUnique <- drugToCellLineDataCgp[drugToCellLineDataCgp$Source.Name %in% 
                                           names(which(table(drugToCellLineDataCgp$Source.Name) == 
                                                         1)), ]
    rownames(pDataUnique) <- pDataUnique$Source.Name
    commonCellLines <- rownames(pDataUnique)[rownames(pDataUnique) %in% 
                                               names(ic50s)]
    pDataUniqueOrd <- pDataUnique[commonCellLines, ]
    ic50sOrd <- ic50s[commonCellLines]
    trainDataOrd <- gdsc_brainarray_syms[, pDataUniqueOrd$Array.Data.File]
    return(list(ic50sOrd = ic50sOrd, trainDataOrd = trainDataOrd))
  }
  else if (dataset == "cgp2016") {
    cat("\nUsing updated CGP 2016 datsets for prediction\n\n")
    if (!tissueType %in% c("all", "aero_digestive_tract", 
                           "blood", "bone", "breast", "digestive_system", "lung", 
                           "nervous_system", "skin", "urogenital_system","allSolidTumors")) 
      stop("ERROR: the tissue type specified must be one of \"all\", \"aero_digestive_tract\", \"blood\", \"bone\", \"breast\", \"digestive_system\", \"lung\", \"nervous_system\", \"skin\", \"urogenital_system\". These tissue types are represented by at least 40 cell lines (although not all 40 may have been screened with each drug).")
    possibleDrugs2016 <- c("Erlotinib", "Rapamycin", "Sunitinib", 
                           "PHA-665752", "MG-132", "Paclitaxel", "Cyclopamine", 
                           "AZ628", "Sorafenib", "VX-680", "Imatinib", "TAE684", 
                           "Crizotinib", "Saracatinib", "S-Trityl-L-cysteine", 
                           "Z-LLNle-CHO", "Dasatinib", "GNF-2", "CGP-60474", 
                           "CGP-082996", "A-770041", "WH-4-023", "WZ-1-84", 
                           "BI-2536", "BMS-536924", "BMS-509744", "CMK", "Pyrimethamine", 
                           "JW-7-52-1", "A-443654", "GW843682X", "MS-275", 
                           "Parthenolide", "KIN001-135", "TGX221", "Bortezomib", 
                           "XMD8-85", "Roscovitine", "Salubrinal", "Lapatinib", 
                           "GSK269962A", "Doxorubicin", "Etoposide", "Gemcitabine", 
                           "Mitomycin C", "Vinorelbine", "NSC-87877", "Bicalutamide", 
                           "QS11", "CP466722", "Midostaurin", "CHIR-99021", 
                           "AP-24534", "AZD6482", "JNK-9L", "PF-562271", "HG-6-64-1", 
                           "JQ1", "JQ12", "DMOG", "FTI-277", "OSU-03012", "Shikonin", 
                           "AKT inhibitor VIII", "Embelin", "FH535", "PAC-1", 
                           "IPA-3", "GSK-650394", "BAY 61-3606", "5-Fluorouracil", 
                           "Thapsigargin", "Obatoclax Mesylate", "BMS-754807", 
                           "Lisitinib", "Bexarotene", "Bleomycin", "LFM-A13", 
                           "GW-2580", "AUY922", "Phenformin", "Bryostatin 1", 
                           "Pazopanib", "LAQ824", "Epothilone B", "GSK1904529A", 
                           "BMS345541", "Tipifarnib", "BMS-708163", "Ruxolitinib", 
                           "AS601245", "Ispinesib Mesylate", "TL-2-105", "AT-7519", 
                           "TAK-715", "BX-912", "ZSTK474", "AS605240", "Genentech Cpd 10", 
                           "GSK1070916", "KIN001-102", "LY317615", "GSK429286A", 
                           "FMK", "QL-XII-47", "CAL-101", "UNC0638", "XL-184", 
                           "WZ3105", "XMD14-99", "AC220", "CP724714", "JW-7-24-1", 
                           "NPK76-II-72-1", "STF-62247", "NG-25", "TL-1-85", 
                           "VX-11e", "FR-180204", "Tubastatin A", "Zibotentan", 
                           "YM155", "NSC-207895", "VNLG/124", "AR-42", "CUDC-101", 
                           "Belinostat", "I-BET-762", "CAY10603", "Linifanib ", 
                           "BIX02189", "CH5424802", "EKB-569", "GSK2126458", 
                           "KIN001-236", "KIN001-244", "KIN001-055", "KIN001-260", 
                           "KIN001-266", "Masitinib", "MP470", "MPS-1-IN-1", 
                           "BHG712", "OSI-930", "OSI-027", "CX-5461", "PHA-793887", 
                           "PI-103", "PIK-93", "SB52334", "TPCA-1", "TG101348", 
                           "Foretinib", "Y-39983", "YM201636", "Tivozanib", 
                           "GSK690693", "SNX-2112", "QL-XI-92", "XMD13-2", 
                           "QL-X-138", "XMD15-27", "T0901317", "EX-527", "THZ-2-49", 
                           "KIN001-270", "THZ-2-102-1", "AICAR", "Camptothecin", 
                           "Vinblastine", "Cisplatin", "Cytarabine", "Docetaxel", 
                           "Methotrexate", "ATRA", "Gefitinib", "Navitoclax", 
                           "Vorinostat", "Nilotinib", "RDEA119", "CI-1040", 
                           "Temsirolimus", "Olaparib", "Veliparib", "Bosutinib", 
                           "Lenalidomide", "Axitinib", "AZD7762", "GW 441756", 
                           "CEP-701", "SB 216763", "17-AAG", "VX-702", "AMG-706", 
                           "KU-55933", "Elesclomol", "Afatinib", "GDC0449", 
                           "PLX4720", "BX-795", "NU-7441", "SL 0101-1", "BIRB 0796", 
                           "JNK Inhibitor VIII", "681640", "Nutlin-3a (-)", 
                           "PD-173074", "ZM-447439", "RO-3306", "MK-2206", 
                           "PD-0332991", "BEZ235", "GDC0941", "AZD8055", "PD-0325901", 
                           "SB590885", "selumetinib", "CCT007093", "EHT 1864", 
                           "Cetuximab", "PF-4708671", "JNJ-26854165", "HG-5-113-01", 
                           "HG-5-88-01", "TW 37", "XMD11-85h", "ZG-10", "XMD8-92", 
                           "QL-VIII-58", "CCT018159", "AG-014699", "SB 505124", 
                           "Tamoxifen", "QL-XII-61", "PFI-1", "IOX2", "YK 4-279", 
                           "(5Z)-7-Oxozeaenol", "piperlongumine", "FK866", 
                           "Talazoparib", "rTRAIL", "UNC1215", "SGC0946", "XAV939", 
                           "Trametinib", "Dabrafenib", "Temozolomide", "Bleomycin (50 uM)", 
                           "SN-38", "MLN4924")
    if (!drug %in% possibleDrugs2016) 
      stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", 
                 (paste(possibleDrugs2016, collapse = ", "))))
    data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
    data(cgp2016ExprRma)
    ic50s <- drugData2016[, "IC50"][drugData2016[, "Drug.name"] == 
                                      drug]
    names(ic50s) <- drugData2016[, "Cell.line.name"][drugData2016[, 
                                                                  "Drug.name"] == drug]
    tissue <- drugData2016[, "Tissue"][drugData2016[, "Drug.name"] == 
                                         drug]
    names(tissue) <- drugData2016[, "Cell.line.name"][drugData2016[, 
                                                                   "Drug.name"] == drug]
    if (tissueType != "all") {
      if (tissueType == "allSolidTumors") {
        tissueType <- ic50s <- ic50s[!(tissue %in% "blood")]
      }
      else {
        ic50s <- ic50s[tissue %in% tissueType]
      }
    }
    commonCellLines <- names(ic50s)[names(ic50s) %in% colnames(cgp2016ExprRma)]
    ic50sOrd <- ic50s[commonCellLines]
    trainDataOrd <- cgp2016ExprRma[, commonCellLines]
    return(list(ic50sOrd = ic50sOrd, trainDataOrd = trainDataOrd))
  }
}

pRRopheticPredict.my <- function (testMatrix, drug, tissueType = "all", batchCorrect = "eb", 
                                  powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2, 
                                  minNumSamples = 10, selection = -1, printOutput = TRUE, 
                                  removeLowVaringGenesFrom = "homogenizeData", dataset = "cgp2014") 
{
  cgpTrainData <- getCGPinfo.my(drug, tissueType, dataset)
  predictedPtype <- calcPhenotype(cgpTrainData$trainDataOrd, 
                                  cgpTrainData$ic50sOrd, testMatrix, batchCorrect = batchCorrect, 
                                  powerTransformPhenotype = powerTransformPhenotype, removeLowVaryingGenes = removeLowVaryingGenes, 
                                  minNumSamples = minNumSamples, selection = selection, 
                                  printOutput = printOutput, removeLowVaringGenesFrom = removeLowVaringGenesFrom)
  return(predictedPtype)
}
