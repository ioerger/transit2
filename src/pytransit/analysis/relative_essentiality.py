import sys

try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True

except Exception as e:
    hasWx = False
    WX_VERSION = 0

if hasWx:
    import wx.xrc
    from wx.lib.buttons import GenBitmapTextButton
    from pubsub import pub
    import wx.adv

import os
import time
import ntpath
import math
import random
import numpy
import scipy.stats
import datetime

from pytransit.analysis import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools
import pytransit.stat_tools as stat_tools



############# GUI ELEMENTS ##################

short_name = "relative_essentiality"
long_name = "Relative Essentiality"
short_desc = "Degree of essentiality of each gene (in a single condition)"
long_desc = """Degree of essentiality of each gene. Ratio (FoldChange) of insertion counts in gene to the global average (1.0=NE, 0.0=ES), with significance (p-value) based on a sampling distribution."""

transposons = ["himar1", "tn5"]
columns = ["ORF","Gene","Annotation","TAsites","Sum","Mean","Sat","NZsites","NZmean","ExpecSum","FoldChange","Pval","Padj"]

class EssentialityAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, EssentialityMethod, EssentialityGUI, [EssentialityFile])




############# FILE ##################

class EssentialityFile(base.TransitFile):

    def __init__(self):
        base.TransitFile.__init__(self, "#RelativeEssentiality", columns)

    def getHeader(self, path):
        DE=0; poslogfc=0; neglogfc=0;
        for line in open(path):
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if float(tmp[-1]) < 0.05:
                DE +=1
                if float(tmp[-3]) > 0:
                    poslogfc+=1
                else:
                    neglogfc+=1

        text = """Results:
    Conditionally - Essentials: %s                     #TRI update this
        Less Essential in Experimental datasets: %s
        More Essential in Experimental datasets: %s
            """ % (DE, poslogfc, neglogfc)
        return text


    def getMenus(self):
        menus = []
        menus.append(("Display in Track View", self.displayInTrackView))
        menus.append(("Display Histogram", self.displayHistogram))
        return menus

    def displayHistogram(self, displayFrame, event):
            gene = displayFrame.grid.GetCellValue(displayFrame.row, 0)
            filepath = os.path.join(ntpath.dirname(displayFrame.path), transit_tools.fetch_name(displayFrame.path))
            filename = os.path.join(filepath, gene+".png")
            if os.path.exists(filename):
                imgWindow = pytransit.fileDisplay.ImgFrame(None, filename)
                imgWindow.Show()
            else:
                transit_tools.ShowError(MSG="Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on.")
                print("Error Displaying File. Histogram image does not exist.")




############# GUI ##################

class EssentialityGUI(base.AnalysisGUI):

    def definePanel(self, wxobj):
        self.wxobj = wxobj
        resamplingPanel = wx.Panel( self.wxobj.optionsWindow, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )

        resamplingSizer = wx.BoxSizer( wx.VERTICAL )

        resamplingLabel = wx.StaticText( resamplingPanel, wx.ID_ANY, u"resampling Options", wx.DefaultPosition, (160,-1), 0 )
        resamplingLabel.SetFont( wx.Font( 10, wx.DEFAULT, wx.NORMAL, wx.BOLD) )
        resamplingSizer.Add( resamplingLabel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        resamplingTopSizer = wx.BoxSizer( wx.HORIZONTAL )

        resamplingTopSizer2 = wx.BoxSizer( wx.HORIZONTAL )

        resamplingLabelSizer = wx.BoxSizer( wx.VERTICAL )

        mainSizer1 = wx.BoxSizer( wx.VERTICAL )

        #(, , Sizer) = self.defineChoiceBox(resamplingPanel, u"", u"", "")
        #mainSizer1.Add(Sizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Samples
        (resamplingSampleLabel, self.wxobj.resamplingSampleText, sampleSizer) = self.defineTextBox(resamplingPanel, u"Samples:", u"10000", "Number of samples to take when estimating the resampling histogram. More samples give more accurate estimates of the p-values at the cost of computation time.")
        mainSizer1.Add(sampleSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Pseudocount
        (resamplingPseudocountLabel, self.wxobj.resamplingPseudocountText, pseudoSizer) = self.defineTextBox(resamplingPanel, u"Pseudocount:", u"0.0", "Adds pseudo-counts to the each data-point. Useful to dampen the effects of small counts which may lead to deceptively high log-FC.")
        mainSizer1.Add(pseudoSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

        # Norm
        resamplingNormChoiceChoices = [ u"TTR", u"nzmean", u"totreads", u'zinfnb', u'quantile', u"betageom", u"nonorm" ]
        (resamplingNormLabel, self.wxobj.resamplingNormChoice, normSizer) = self.defineChoiceBox(resamplingPanel, u"Normalization: ", resamplingNormChoiceChoices, "Choice of normalization method. The default choice, 'TTR', normalizes datasets to have the same expected count (while not being sensative to outliers). Read documentation for a description other methods. ")
        mainSizer1.Add(normSizer, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )




        resamplingSizer.Add( mainSizer1, 1, wx.EXPAND, 5 )






        # LOESS Check
        (self.wxobj.resamplingLoessCheck, loessCheckSizer) = self.defineCheckBox(resamplingPanel, labelText="Correct for Genome Positional Bias", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Check to correct read-counts for possible regional biase using LOESS. Clicking on the button below will plot a preview, which is helpful to visualize the possible bias in the counts.")
        resamplingSizer.Add( loessCheckSizer, 0, wx.EXPAND, 5 )

        # LOESS Button
        self.wxobj.resamplingLoessPrev = wx.Button( resamplingPanel, wx.ID_ANY, u"Preview LOESS fit", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingSizer.Add( self.wxobj.resamplingLoessPrev, 0, wx.ALL|wx.CENTER, 5 )

        # Adaptive Check
        (self.wxobj.resamplingAdaptiveCheckBox, adaptiveSizer) = self.defineCheckBox(resamplingPanel, labelText="Adaptive Resampling (Faster)", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Dynamically stops permutations early if it is unlikely the ORF will be significant given the results so far. Improves performance, though p-value calculations for genes that are not differentially essential will be less accurate.")
        resamplingSizer.Add( adaptiveSizer, 0, wx.EXPAND, 5 )

        # Histogram Check
        (self.wxobj.resamplingHistogramCheckBox, histSizer) = self.defineCheckBox(resamplingPanel, labelText="Generate Resampling Histograms", widgetCheck=False, widgetSize=(-1,-1), tooltipText="Creates .png images with the resampling histogram for each of the ORFs. Histogram images are created in a folder with the same name as the output file.")
        resamplingSizer.Add(histSizer, 0, wx.EXPAND, 5 )


        # Zeros Check
        (self.wxobj.resamplingZeroCheckBox, zeroSizer) = self.defineCheckBox(resamplingPanel, labelText="Include sites with all zeros", widgetCheck=True, widgetSize=(-1,-1), tooltipText="Includes sites that are empty (zero) across all datasets. Unchecking this may be useful for tn5 datasets, where all nucleotides are possible insertion sites and will have a large number of empty sites (significantly slowing down computation and affecting estimates).")
        resamplingSizer.Add(zeroSizer, 0, wx.EXPAND, 5 )


        resamplingButton = wx.Button( resamplingPanel, wx.ID_ANY, u"Run resampling", wx.DefaultPosition, wx.DefaultSize, 0 )
        resamplingSizer.Add( resamplingButton, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


        resamplingPanel.SetSizer( resamplingSizer )
        resamplingPanel.Layout()
        resamplingSizer.Fit( resamplingPanel )

        #Connect events
        resamplingButton.Bind( wx.EVT_BUTTON, self.wxobj.RunMethod )
        self.wxobj.resamplingLoessPrev.Bind(wx.EVT_BUTTON, self.wxobj.LoessPrevFunc)

        self.panel = resamplingPanel

    def GlobalEnable(self):
        self.wxobj.ctrlLibText.Enable()
        self.wxobj.expLibText.Enable()

    def GlobalDisable(self):
        self.wxobj.ctrlLibText.Disable()
        self.wxobj.expLibText.Disable()



########## CLASS #######################

class EssentialityMethod(base.DualConditionMethod):
    """
    resampling

    """
    def __init__(self,
                ctrldata,
                #expdata,        #TRI
                annotation_path,
                output_file,
                normalization="TTR",
                samples=10000,
                adaptive=False,
                doHistogram=False,
                includeZeros=False,
                pseudocount=0.0,
                replicates="Sum",
                LOESS=False,
                ignoreCodon=True,
                NTerminus=0.0,
                CTerminus=0.0,
                ctrl_lib_str="",
                exp_lib_str="",
                wxobj=None, Z = False, diffStrains = False, annotation_path_exp = "", combinedWigParams = None):

        base.DualConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, ctrldata , None , annotation_path, output_file, normalization=normalization, replicates=replicates, LOESS=LOESS, NTerminus=NTerminus, CTerminus=CTerminus, wxobj=wxobj)

        self.Z = Z
        self.samples = samples
        self.adaptive = adaptive
        self.doHistogram = doHistogram
        self.includeZeros = includeZeros
        self.pseudocount = pseudocount
        self.ctrl_lib_str = ctrl_lib_str
        self.exp_lib_str = exp_lib_str
        self.diffStrains = diffStrains
        self.annotation_path_exp = annotation_path_exp if diffStrains else annotation_path
        self.combinedWigParams = combinedWigParams

    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get Annotation file
        annot_paths = wxobj.annotation.split(",")
        annotationPath = annot_paths[0]
        diffStrains = False
        annotationPathExp = ""
        if len(annot_paths) == 2:
            annotationPathExp = annot_paths[1]
            diffStrains = True

        if not transit_tools.validate_annotation(annotationPath):
            return None

        if annotationPathExp and not transit_tools.validate_annotation(annotationPathExp):
            return None

        #Get selected files
        ctrldata = wxobj.ctrlSelected()
        expdata = wxobj.expSelected()
        if not transit_tools.validate_both_datasets(ctrldata, expdata):
            return None

        #Validate transposon types
        if not transit_tools.validate_transposons_used(ctrldata+expdata, transposons):
            return None


        #Read the parameters from the wxPython widgets
        ignoreCodon = True
        samples = int(wxobj.resamplingSampleText.GetValue())
        normalization = wxobj.resamplingNormChoice.GetString(wxobj.resamplingNormChoice.GetCurrentSelection())
        replicates="Sum"
        adaptive = wxobj.resamplingAdaptiveCheckBox.GetValue()
        doHistogram = wxobj.resamplingHistogramCheckBox.GetValue()

        includeZeros = wxobj.resamplingZeroCheckBox.GetValue()
        pseudocount = float(wxobj.resamplingPseudocountText.GetValue())
        LOESS = wxobj.resamplingLoessCheck.GetValue()

        # Global Parameters
        NTerminus = float(wxobj.globalNTerminusText.GetValue())
        CTerminus = float(wxobj.globalCTerminusText.GetValue())
        ctrl_lib_str = wxobj.ctrlLibText.GetValue()
        exp_lib_str = wxobj.expLibText.GetValue()


        #Get output path
        defaultFileName = "resampling_output_s%d_pc%1.2f" % (samples, pseudocount)
        if adaptive: defaultFileName+= "_adaptive"
        if includeZeros: defaultFileName+= "_iz"
        defaultFileName+=".dat"

        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")


        return self(ctrldata,
                expdata,
                annotationPath,
                output_file,
                normalization,
                samples,
                adaptive,
                doHistogram,
                includeZeros,
                pseudocount,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus,
                ctrl_lib_str,
                exp_lib_str, wxobj, Z = False, diffStrains = diffStrains, annotation_path_exp = annotationPathExp)

    @classmethod
    def fromargs(self, rawargs):

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        isCombinedWig = True if kwargs.get('c', False) else False
        combinedWigParams = None
        if isCombinedWig:
            if (len(args) != 5):
                print("Error: Incorrect number of args. See usage")
                print(self.usage_string())
                sys.exit(0)
            combinedWigParams = {
                "combined_wig": kwargs.get('c'),
                "samples_metadata": args[0],
                "conditions": [args[1].lower(), args[2].lower()]
            }
            annot_paths = args[3].split(",")
            ctrldata = ""
            expdata = ""
            output_path = args[4]
        else:
            if (len(args) != 3):
                print("Error: Incorrect number of args. See usage")
                print(self.usage_string())
                sys.exit(0)
            ctrldata = args[0].split(",")
            #expdata = args[1].split(",")
            annot_paths = args[1].split(",") #TRI
            output_path = args[2] #TRI
        annotationPath = annot_paths[0]
        diffStrains = False
        annotationPathExp = ""
        if len(annot_paths) == 2:
            annotationPathExp = annot_paths[1]
            diffStrains = True
        if (diffStrains and isCombinedWig):
            print("Error: Cannot have combined wig and different annotation files.")
            sys.exit(0)

        output_file = open(output_path, "w")

        normalization = kwargs.get("n", "TTR")
        samples = int(kwargs.get("s", 10000))
        adaptive = kwargs.get("a", False) #TRI eventually can get rid of these
        doHistogram = kwargs.get("h", False)
        replicates = kwargs.get("r", "Sum")
        excludeZeros = kwargs.get("ez", False)
        includeZeros = not excludeZeros
        pseudocount = float(kwargs.get("pc", 0.00))

        Z = True if "Z" in kwargs else False

        LOESS = kwargs.get("l", False)
        ignoreCodon = True

        NTerminus = float(kwargs.get("iN", 0.00))
        CTerminus = float(kwargs.get("iC", 0.00))
        ctrl_lib_str = kwargs.get("-ctrl_lib", "")
        exp_lib_str = kwargs.get("-exp_lib", "")

        return self(ctrldata,
                #expdata,
                annotationPath,
                output_file,
                normalization,
                samples,
                adaptive,
                doHistogram,
                includeZeros,
                pseudocount,
                replicates,
                LOESS,
                ignoreCodon,
                NTerminus,
                CTerminus,
                ctrl_lib_str,
                exp_lib_str, 
                Z = Z, diffStrains = diffStrains, annotation_path_exp = annotationPathExp, combinedWigParams = combinedWigParams)

    def preprocess_data(self, position, data):
        (K,N) = data.shape

        if self.normalization != "nonorm":
            self.transit_message("Normalizing using: %s" % self.normalization)
            #(data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata+self.expdata, self.annotation_path)
            (data, factors) = norm_tools.normalize_data(data, self.normalization, self.ctrldata, self.annotation_path)

        if self.LOESS:
            self.transit_message("Performing LOESS Correction")
            for j in range(K):
                data[j] = stat_tools.loess_correction(position, data[j])

        return data

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, None) for f in filenamesInCombWig]

    def filter_wigs_by_conditions(self, data, conditions, included_conditions):
        """
            Filters conditions from wig to ctrl, exp conditions only
            ([[Wigdata]], [ConditionCtrl, ConditionExp]) -> Tuple([[Wigdata]], [Condition])
        """
        d_filtered, cond_filtered = [], []
        if len(included_conditions) != 2:
            self.transit_error("Only 2 conditions expected", included_conditions)
            sys.exit(0)
        for i, c in enumerate(conditions):
          if (c.lower() in included_conditions):
            d_filtered.append(data[i])
            cond_filtered.append(conditions[i])

        return (numpy.array(d_filtered), numpy.array(cond_filtered))

    # remove all runs of zeros of length >= W

    def remove_essential_regions(self, wig,W):
      runs = []
      i,n = 0,len(wig)
      while i<n:
        if wig[i]>0: i += 1
        else:
          j = i
          while j<n and wig[j]==0: j += 1
          runs.append((i,j))
          i = j
      counts = []
      for k,(i,j) in enumerate(runs):
        if j-i<W: counts += wig[i:j]
        if k<len(runs)-1:
          next = runs[k+1][0]
          counts += wig[j:next]
      return counts

    # given k pools of counts, draw n from each and concatenate them; repeat N times

    def sample_count_pools(self,count_pools,size,times):
      samples = []
      for i in range(times):
        sample = []
        for counts in count_pools: sample += random.sample(counts,size) # without replacement
        samples.append(sample) 
      return samples

    def Run(self):

        #if not self.wxobj:
        #    # Force matplotlib to use good backend for png.
        #    import matplotlib.pyplot as plt
        #elif "matplotlib.pyplot" not in sys.modules:
        try:
            import matplotlib.pyplot as plt
        except:
            print("Error: cannot do histograms")
            self.doHistogram = False


        self.transit_message("Starting Relative Essentiality Method")
        start_time = time.time()

        #Get orf data
        self.transit_message("Getting Data")
        if self.diffStrains:
            self.transit_message("Multiple annotation files found")
            self.transit_message("Mapping ctrl data to {0}, exp data to {1}".format(self.annotation_path, self.annotation_path_exp))

        if self.combinedWigParams:
            (position, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combinedWigParams['combined_wig'])
            conditionsByFile, _, _, _ = tnseq_tools.read_samples_metadata(self.combinedWigParams['samples_metadata'])
            conditions = self.wigs_to_conditions(conditionsByFile, filenamesInCombWig)
            data, conditions = self.filter_wigs_by_conditions(data, conditions, self.combinedWigParams['conditions'])
            data_ctrl = numpy.array([d for i, d in enumerate(data) if conditions[i].lower() == self.combinedWigParams['conditions'][0]])
            position_ctrl, position_exp = position, position
        else:
            (data_ctrl, position_ctrl) = transit_tools.get_validated_data(self.ctrldata, wxobj=self.wxobj)
        (K_ctrl, N_ctrl) = data_ctrl.shape

#TRI handle diff num of sites
#        if not self.diffStrains and (N_ctrl != N_exp):

        self.transit_message("Preprocessing Ctrl data...")
        data_ctrl = self.preprocess_data(position_ctrl, data_ctrl) # call normalization

        G_ctrl = tnseq_tools.Genes(self.ctrldata, self.annotation_path, ignoreCodon=self.ignoreCodon, nterm=self.NTerminus, cterm=self.CTerminus, data=data_ctrl, position=position_ctrl)

        self.transit_message("creating pools of counts in non-essential regions")
        count_pools = []
        for i,counts in enumerate(data_ctrl):
          counts = list(counts)
          n,sat = len(counts),tnseq_tools.saturation(counts)
          R = tnseq_tools.ExpectedRuns(n,1.0-sat) # pnon
          varR = tnseq_tools.VarR(n,1.0-sat) # pnon
          W = R-varR
          pool = self.remove_essential_regions(counts,W)
          print("%s: saturation=%0.1f%%, E[R]=%0.1f, var(R)=%0.1f, run-length for filtering out essential regions=%0.1f, pool: sites=%s, sat=%0.1f%%" % (self.ctrldata[i],100.*sat,R,varR,W,len(pool),100.*tnseq_tools.saturation(pool)))
          count_pools.append(pool)

        results = self.run_essentiality(G_ctrl,count_pools)
        self.write_output(results, start_time)

        self.finish()
        self.transit_message("Finished Relative Essentiality Method")

    def write_output(self, results, start_time):

        self.output.write("#RelativeEssentiality\n")
        if self.wxobj:
            members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
            memberstr = ""
            for m in members:
                memberstr += "%s = %s, " % (m, getattr(self, m))
            self.output.write("#GUI with: norm=%s, samples=%s, pseudocounts=%1.2f, adaptive=%s, histogram=%s, includeZeros=%s, output=%s\n" % (self.normalization, self.samples, self.pseudocount, self.adaptive, self.doHistogram, self.includeZeros, self.output.name.encode('utf-8')))
        else:
            self.output.write("#Console: python %s\n" % " ".join(sys.argv))
        self.output.write("#Control Data: %s\n" % (",".join(self.ctrldata).encode('utf-8')))
        self.output.write("#Annotation path: %s %s\n" % (self.annotation_path.encode('utf-8'), self.annotation_path_exp.encode('utf-8') if self.diffStrains else ''))
        self.output.write("#Time: %s\n" % (time.time() - start_time))

        global columns # consider redefining columns above (for GUI)
        self.output.write("#%s\n" % "\t".join(columns))

        for row in results: self.output.write('\t'.join([str(x) for x in row])+"\n")
        self.output.close()

        self.transit_message("Adding File: %s" % (self.output.name))
        self.add_file(filetype="RelativeEssentiality")

    def run_essentiality(self, G_ctrl, count_pools):
        results = []
        Ngenes = len(G_ctrl)
        count = 0
        self.progress_range(Ngenes)
        cache = {}

        for gene in G_ctrl:
            count+=1
            if gene.n==0: continue
            cnts = gene.reads.flatten() #TRI consider adding self.pseudocounts
            n = gene.n

            nonzeros = [cnts[x] for x in numpy.nonzero(cnts)[0]]
            zeros,NZmean = n-len(nonzeros),numpy.mean(nonzeros) if len(nonzeros)>0 else 0
            sat = len(nonzeros)/float(n)
            tot = sum(cnts)
            mn = tot/float(n)

            N,alpha = 10000,0.05
            if n in cache: sample = cache[n]
            else: 
              sample = self.sample_count_pools(count_pools,n,N) # returns a list of N lists of length k*n, where K in num of pools, and n is num of TA sites
              cache[n] = sample
            samplesums = [sum(lst) for lst in sample]
            meansum = numpy.mean(samplesums)
            PC = 1
            rel = (sum(cnts)+PC)/float(meansum+PC)
            LFC = numpy.log2(rel)

            lesser = len(list(filter(lambda x: x<=tot,samplesums)))
            greater = len(list(filter(lambda x: x>=tot,samplesums)))
            pval = min(lesser,greater)/float(N)

            vals = [gene.orf, gene.name, gene.desc, gene.n]
            vals += ["%s" % int(tot),"%0.1f" % mn,"%0.3f" % sat,len(nonzeros),"%0.1f" % NZmean]
            vals += [int(meansum),"%0.3f" % rel,pval]

            #print('\t'.join([str(x) for x in vals]))
            results.append(vals)
            
            # Update progress
            text = "Running Relative Essentiality Method... %5.1f%%" % (100.0*count/Ngenes)
            self.progress_update(text, count)


        self.transit_message("") # Printing empty line to flush stdout
        self.transit_message("Performing Benjamini-Hochberg Correction")
        results.sort()
        qvals = stat_tools.BH_fdr_correction([row[-1] for row in results])
        results = [x+[y] for x,y in zip(results,qvals)]
        return results

    @classmethod
    def usage_string(self):
        return """
        python %s relative_essentiality <comma-separated-list-of-wig-files> <annotation .prot_table or GFF3> <output file> [Optional Arguments]

        Optional Arguments:
        -s <integer>    :=  Number of samples. Default: -s 10000
        -n <string>     :=  Normalization method. Default: -n TTR
        -pc             :=  Pseudocounts to be added at each site.
        -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 0
        -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 0

        """ % sys.argv[0]

if __name__ == "__main__":

    (args, kwargs) = transit_tools.cleanargs(sys.argv)

    #TODO: Figure out issue with inputs (transit requires initial method name, running as script does not !!!!)
    
    G = ResamplingMethod.fromargs(sys.argv[1:])

    G.console_message("Printing the member variables:")
    G.print_members()

    print("")
    print("Running:")

    G.Run()


