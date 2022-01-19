mod_Presentation_ui <- function(id) {
  ns <- NS(id)
  # Presentation

  tabItem(tabName = "Presentation",
  
    fluidPage(
      h1("Flinovo project"),
      
      fluidRow(
        column(8,
          h2("Introduction"),
          
          p("The objective of the FLINOVO project is to analyze the clinical response to treatment of patients with follicular lymphoma and to draw correlations between responder and non-responder groups of patients responders to RCHOP.
                The generation of bioinformatics data related to the analyses of interest depends upstream on the implementation of a biological experimental protocol subject to the use of a particular animal model model: the AVI-PDX model. 
                This model, designed by the company Oncofactory, consists in the very rapid creation of miniaturized replicates of primary tumors for each patient, made by xenogree of tumor cells in chicken embryos. 
                The therapeutic molecules (RCHOP or excipient) are delivered to the embryo for direct evaluation of their efficacy. This technique allows to mimic the heterogeneity of patient tumors, their responses to therapies and 
                their interactions with the microenvironment. Once the cells have been transplanted and subjected to treatment for 24 hours, an enzymatic dissociation step is performed to extract the cells of interest and analyze their status."),
          
          HTML("Thus, for each patient we analyze 3 different cell populations: 
                <ul>
                  <li>the control condition: a first one not subjected to grafting</li>
                  <li>the excipient condition: a second one subjected to grafting but not to treatment</li>
                  <li>the RCHOP condition: a last one subjected to both grafting and treatment</li>
                </ul>
                This experimental protocol is applied in order to distinguish the effects specific to the graft and specific to the treatment. In order to analyze in a single experiment the three experimental conditions, 
                the cells will be multiplexed, either independently labeled with a barcoded antibody carrying a molecular tag, then sorted by flow cytometry to generate libraries.<br><br>"),

          
          
          h2("Bioinformatics Workflow"),
          p(""),
          
          
          h3("Patient by patient."),
          
          h4("I. Pre-processing."),
          
          HTML("<li>Exctraction of the transcriptomic information : Long story short.</li>
                    <ul>
                      <li>1. CellRanger mkfastq : Convertion from BCL raw sequencing files to FASTQ files</li>
                      <li>2. CellRanger Count : mRNA library from GRCh38 human genome references</li>
                      <li>3. CITE-SEQ-count : To demultiplex every condition merged inside our sequenced sample </li>
                      <li>4. CellRanger VDJ : To analyse every VDJ fragment of every sample </li>
                    </ul>"),
          
          HTML("<li>Construction of the <a href='https://satijalab.org/seurat/'>Seurat</a> object :</li>
                    <ul>
                      <li>Construction of the mRNA array from CellRangerCount output as a seurat object</li>
                      <li>Construction of the HTO array from CITE-SEQ-count output as a seurat object</li>
                      <li>Merge of mRNA and HTO arrays to form a final seurat object, and first metadata cleaning</li>
                      <li>Addition of percentages for mitochondrial (MT), ribosomal (RP[SL]) and immunoglobulin (IG) gene expression</li>
                    </ul>"),
          
          HTML("<li>Normalisation :</li>
                    <ul>
                      <li>After many test, we choose to apply the SCTransform normalisation method from this <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1'>paper</a>.
                      More information on this <a href='https://github.com/satijalab/sctransform'>Github</a> repository.</li>
                      <li>Applied method : glmGamPoi</li>
                    </ul>"),
          
          HTML("<li>Way of visualisation added :</li>
                    <ul>
                      <li>PCA : standards options</li>
                      <li>UMAP : 40 first dimensions, with the 'umap-learn' method</li>
                      <li>FindNeighbors : 40 first dimensions</li>
                      <li>FindClusters : with 0.6 as resolution</li>
                      <li>TSNE : 40 first dimensions, with 100 as perplexity</li>
                    </ul>"),
          
          br(),
          h4("II. Addition of metadata."),
          HTML("<li>Cell phenotype identification :</li>
                      <ul>
                        <li>To perform the cell phenotype identification, we used the MonacoImmuneData database from the <a href='https://bioconductor.org/packages/devel/data/experiment/manuals/celldex/man/celldex.pdf'>celldex</a> R package.</li>
                        <li>Both 'main' and 'fine' labels were added as metadata to characterize the cell phenotype.</li>
                      </ul>
                    <li>Analyses of BCR receptors in polyclonal B cell populations : </li></ul>
                      <ul>
                        <li>Data extracted from matrixs and BCL files of the CellRanger VDJ output</li>
                        <li>Visualisation of the VDJC genes expressions, of each expressed chain (heavy and light) and of the CDR3.</li>
                      </ul>
                    <li>Genotype of Transcription : GoT</li></ul>
                      <ul>
                        <li>To observe the presence or not of 10 different mutations shared between 5 hotspots and distributed on 3 genes :</li>
                        <ul>
                          <li>BCL2 : K22K, L23L</li>
                          <li>CD79B : Y196H</li>
                          <li>EZH2 : Y146C, Y146F, Y146H, Y146N, Y146S, A192G, A192V</li>
                        </ul>
                        <li>Elaboration of a RNA librariy from read spcificly expressed by those area of genes : 100 UMIs without replicates was was extracted via python for each cell .</li>
                        <li>Catégorisation by the <a href='https://github.com/landau-lab/IronThrone-GoT'>Irothrone-GoT software</a> of the cell's genotypes as MUT, WT or AMB.</li>
                        <li>A cell is considered to be WT or MUT if and only if more than 80% of its IMUs support the same genotype. below this threshold the cells are scored as AMB.</li>
                      </ul>
                    <li>Phase cycle information : CellCycleScoring function form the Seurat package</li></ul>"),
          
          br(),
          h4("III. Additional analyses."),
          
          HTML("<li>Velocity :</li>
                    <ul>
                      <li>The following python software <a href='http://velocyto.org/velocyto.py/index.html'>velocyto run10x</a> were runned directly from the output of CellRanger Count to determine the velocity of each sample.</li>
                      <li>The .loom file output is converted as a seurat object.</li>
                      <li>Two ways of calcul were programed to show the cell velocity : </li>
                        <ul><li>'Seurat independant' calcul : Datas are normalised by SCTransform from the 'spliced' assay and the function RunVelocity is applied. A new seurat object is built.</li></ul>
                        <ul><li>'Seurat dependant' calcul : The assays spliced/unspliced/ambiguous are merged with the RNA/HTO/SCT assay from the previous initial seurat object and the function RunVelocity is applied, without 'SCT' renormalisation.</li></ul>
                      <li>The restuls are shown by the show.velocity.on.embedding.cor function from the <a href='https://github.com/velocyto-team/velocyto.R'>velocito.R</a> package.</li>
                   </ul>"),
          
          HTML("<li>Entropy :</li>
                    <ul><li>Intercellular entropy : endowing a gene with an entropy value</li>
                      <ul>
                        <li>Bootstrap's utilisation : 100 random counts's cells with replace are taken</li>
                        <li>The entropy is calculated for each gene of the matrix by the estimateur_bub_entropie_discrete function</li>
                        <li>Every entropy is determined at least 50 time and the mean entropy is calculated for each gene</li>
                        <li>Boxplot : RCHOP / Excipient</li>
                      </ul>
                      <li>Intracellular : endowing a cell with an entropy value</li>
                    </ul>"),
          
          HTML("<li><a href='https://github.com/vitkl/ParetoTI'>Pareto</a> : Archetypal analysis and task inference </li>
                    <ul>
                      <li>This analyse is based on this <a href='https://vitkl.github.io/ParetoTI/articles/Hepatocyte_example.html#load-data-from-geo-and-filter-as-described-in-the-paper-normalise-and-pcs-for-finding-polytopes'>tuto</a>. See more information <a href='https://www.weizmann.ac.il/mcb/UriAlon/download/ParTI'>here</a>.</li>
                      <li>The initial seurat object is converted to a SCE object</li>
                      <li>Load data from GEO : Annotation MT genes</li>
                      <li>Filtration as described in the paper : Cells with more less than 1000 or more than 30000 UMI, genes and cells with too many zeros are filtered</li>
                      <li>Then we normalise the gene expression by cell sum factors and log-transformation is applied</li>
                      <li>Find archetypes : Examination of the polytope with best k & look at known markers of subpopulations</li>
                      <li>Keep the N archetypes which's explained more than 80% of the variability</li>
                      <li>Find genes and gene sets enriched near vertices</li>
                      <li>Map GO annotations and measure activities according to the archetype : Cellular Composant + Biolobical Process + Molecular Function</li>
                      <li>Take a look at top genes and functions for each archetype</li>
                   </ul>"),
          #<a href=''></a> 
          
          br(),
          h4("IV. Data exploration."),
          HTML("<ul>
                      <li>Subset of dataset to split combinaison of post-greffe, pre-greffe, B & T cells or just B cells.</li>
                      <li>Application of a diet function to deal with light dataset inside this application.</li>
                   </ul>"),
          
          
          br(),
          h3("Meta analyses : Integration of all the patients."),
          HTML("<li>Construction of Seurat objects for each patient as discribed before and merging every one of them to a merged final seurat object</li>"),
          HTML("<li>Addition the phenotypical information as discribed before</li>"),
          HTML("<li>Integration : </li>
                          <ul>
                            <li>Split the merged object by patient to a list</li>
                            <li>Normalisation by SCTranform</li>
                            <li>A first set of gene are selected as 'integration features genes' : refect the differents feature between every datasets</li>
                            <li>A second set of gene are selected as 'anchors features genes' : refect the common points between every datasets </li>
                            <li>Integration by seurat</li>
                         </ul>"),
          
          HTML("<li>Way of visualisation added : As discribed before but with every way of visualisation focused on the variable features genes set calculted from the SCT assays without every gene related to the immunoglobuline gene expression.</li>"),
          HTML("<li>Addition of metadatas</li>"),
          
          HTML("<li>Differential gene expression : With Seurat</li>
                          <ul>
                            <li>Every genes from the RNA data matrix are taken to calculate gene unless every one of them related to the immunoglobuline gene expression. (regex : IG[HKL][VDJ]|IGHG[1-4]|IGH[MDE]|IGKC|IGLL|IGLC[1-7]|IGHA[1-2]|TR[ABGD][CV]) </li>
                            <li>Function FindMarkers, assay : 'RNA', slot : 'data', probability law : negbinom </li>
                          </ul>"),
          
          HTML("This application has been developed with the R software, using the <a href='http://shiny.rstudio.com'>Shiny</a> package."),

        ),
        column(4,
          box(title = "Sunburst : Metadata", plotlyOutput(ns("sunburst"), width = "100%",  height = "800px"), collapsible = T, width = 12 ),
          box(title = "Histogram summary", plotOutput(ns("hist"), width = "100%",  height = "800px"), collapsible = T, width = 12 )
          #box(title = "PCA", plotOutput(ns("PCA"), width = "100%",  height = "350px"), collapsible = T, width = 12 ),
          #box(title = "UMAP", plotOutput(ns("UMAP"), width = "100%",  height = "350px"), collapsible = T, width = 12 ),
          #box(title = "TSNE", plotOutput(ns("TSNE"), width = "100%",  height = "350px"), collapsible = T, width = 12 )
        )
      )
    )
  )
}

mod_Presentation_server <- function(input, output, session, r) {
  ns <- session$ns
  
  load(file="inst/app/www/presentation.RData")
  load(file="inst/app/www/presentation.RData")
  
  output$sunburst <- renderPlotly({
    plot_ly(sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3")), type = 'sunburst', #FDB462
          branchvalues = 'total', hoverinfo = "text", hovertext = paste(sunburst$labels, ":", round((sunburst$values/sunburst$values[1])*100,2),"%", "\nTotal : " , sunburst$values))
  })

  # output$PCA <- renderPlot({Seurat::DimPlot(object = all, label.size = 5, pt.size = 0.3, reduction = 'pca', label = TRUE) & Seurat::NoLegend() & Seurat::NoAxes()})
  # output$UMAP <- renderPlot({Seurat::DimPlot(object = all, label.size = 5, pt.size = 0.3, reduction = 'umap', label = TRUE) & Seurat::NoLegend() & Seurat::NoAxes()})
  # output$TSNE <- renderPlot({Seurat::DimPlot(object = all, label.size = 5, pt.size = 0.3, reduction = 'tsne', label = TRUE) & Seurat::NoLegend() & Seurat::NoAxes()})
  
  output$hist <- renderPlot({ggplot(data=result, aes(x=Sample, y=Number, fill=Condition)) + geom_bar(stat="identity", position=position_dodge())+scale_fill_brewer(palette="Paired") + theme_minimal()})
  
}

## To be copied in the UI
# mod_Presentation_ui("Presentation_ui_1")

## To be copied in the server
# callModule(mod_Presentation_server, "Presentation_ui_1")