# ============================================================================
# 00_theme_colors.R
# Project: Decoding Immune Landscapes and Therapeutic Targets in Intracerebral Hemorrhage
# Author:  Dimitrios Kyriakis
# ============================================================================

# ==============================================================================
# Shared color palettes used across Figure 1 and Figure 2 scripts.
# Source after 00_setup.R.
# ==============================================================================

# ---- Cell population / subtype colors (Figure 1: Immune CellType; Figure 2 labels) ----
Project_Colors <- c(
  'Astrocytes' = '#D9D9D9',
  'B Cells'    = '#FDB462',
  'Myeloid'    = '#FB8072',
  'MONO'       = '#80B1D3',
  'Murel'      = '#BEBADA',
  'NEUT'       = '#B3DE69',
  'Oligo'      = '#FCCDE5',
  'T Cells'    = '#8DD3C7',
  'MONO/NEUT'  = '#80B1D3',

  'MG_Homeo' = '#279e68',
  'MG_PVM'   = '#ff7f0e',
  'MG_Prolif'= '#d62728',
  'MG_Adapt' = '#1f77b4',
  'MG_exAM'  = '#ff9896',

  'MG_Homeo_FRMD4A' = '#279e68',
  'MG_Homeo_PICALM' = '#98df8a',

  'MG_PVM_CD163' = '#ff7f0e',
  'MG_PVM_GPNMB' = '#ffbb78',

  'MG_Prolif_MKI67' = '#d62728',

  'MG_Adapt_TMEM163' = '#1f77b4',
  'MG_Adapt_AIF1'    = '#c5b0d5',
  'MG_Adapt_IFI44L'  = '#17becf',
  'MG_Adapt_HIF1A'   = '#aec7e8',
  'MG_Adapt_CCL3'    = '#aa40fc',
  'MG_Adapt_HSPA1A'  = '#9edae5',
  'MG_Adapt_HIST'    = '#dbdb8d',

  'MG_exAM_ERN1' = '#ff9896'
)

# ---- Figure 2: MTC subclass colors (7 metacell clusters) ----
MTC_Colors <- c(
  "MTC_1"      = "#D55E00",
  "MTC_2"      = "#E69F00",
  "MTC_3"      = "#A6611A",
  "MTC_4"      = "#92C5DE",
  "MTC_5"      = "#053061",
  "MTC_6"      = "#2166AC",
  "MTC_Prolif" = "#dbdb8d"
)

# ---- Figure 2: MTC_123 vs MTC_456 two-class colors (from Set3/RdGy palettes) ----
colors_set3   <- brewer.pal(n = 9, name = "Set3")
colors_set3[2] <- colors_set3[8]
colors_set3[8] <- "azure3"
colors_rb     <- rev(brewer.pal(n = 6, name = "RdGy"))
TwoClass_Colors <- c("MTC-456" = colors_rb[6], "MTC-123" = colors_rb[2])

# ---- Figure 2 (supplementary): IREA cytokine family colors ----
Family_Colors <- c(
  'Interferon'                       = '#20aa8b',
  'IL1'                              = '#922a6e',
  'Common gamma chain / IL-13/TSLP'  = '#F6C76D',
  'Common beta chain'                = '#8668c0',
  'IL-6 / IL-12'                     = '#198591',
  'IL10'                             = '#df70b0',
  'IL17'                             = '#096954',
  'Growth Factor'                    = '#472d8c',
  'TNF'                              = '#53b8c4',
  'Complement'                       = '#b7572e',
  'Other'                            = '#91a160'
)

Cytokine_Family <- c(
  'IFN-α1'='Interferon','IFN-β'='Interferon','IFN-ε'='Interferon',
  'IFN-κ'='Interferon','IFN-γ'='Interferon','IFN-λ2'='Interferon',

  'IL-1α'='IL1','IL-1β'='IL1','IL-1Ra'='IL1','IL-18'='IL1',
  'IL-33'='IL1','IL-36α'='IL1','IL-36Ra'='IL1',

  'IL-2'='Common gamma chain / IL-13/TSLP','IL-4'='Common gamma chain / IL-13/TSLP',
  'IL-13'='Common gamma chain / IL-13/TSLP','IL-15'='Common gamma chain / IL-13/TSLP',
  'IL-7'='Common gamma chain / IL-13/TSLP','TSLP'='Common gamma chain / IL-13/TSLP',
  'IL-9'='Common gamma chain / IL-13/TSLP','IL-21'='Common gamma chain / IL-13/TSLP',

  'IL-3'='Common beta chain','IL-5'='Common beta chain','GM-CSF'='Common beta chain',

  'IL-6'='IL-6 / IL-12','IL-11'='IL-6 / IL-12','IL-27'='IL-6 / IL-12',
  'IL-30'='IL-6 / IL-12','IL-31'='IL-6 / IL-12','LIF'='IL-6 / IL-12',
  'OSM'='IL-6 / IL-12','CT-1'='IL-6 / IL-12','NP'='IL-6 / IL-12',
  'IL-12'='IL-6 / IL-12','IL-23'='IL-6 / IL-12','IL-Y'='IL-6 / IL-12',

  'IL-10'='IL10','IL-19'='IL10','IL-20'='IL10','IL-22'='IL10','IL-24'='IL10',

  'IL-17A'='IL17','IL-17B'='IL17','IL-17C'='IL17',
  'IL-17D'='IL17','IL-17E'='IL17','IL-17F'='IL17',

  'FLT3L'='Growth Factor','IL-34'='Growth Factor','M-CSF'='Growth Factor',
  'G-CSF'='Growth Factor','SCF'='Growth Factor','EGF'='Growth Factor',
  'VEGF'='Growth Factor','FGF-β'='Growth Factor','HGF'='Growth Factor','IGF-1'='Growth Factor',

  'LT-α1/β2'='TNF','LT-α2/β1'='TNF','TNF-α'='TNF','OX40L'='TNF',
  'CD40L'='TNF','FasL'='TNF','CD27L'='TNF','CD30L'='TNF','4-1BBL'='TNF',
  'TRAIL'='TNF','RANKL'='TNF','TWEAK'='TNF','APRIL'='TNF','BAFF'='TNF',
  'LIGHT'='TNF','TL1A'='TNF','GITRL'='TNF',

  'C3a'='Complement','C5a'='Complement',

  'PRL'='Other','Leptin'='Other','AdipoQ'='Other','ADSF'='Other',
  'TGF-β1'='Other','GDNF'='Other','PSPN'='Other','Noggin'='Other',
  'Decorin'='Other','TPO'='Other'
)
