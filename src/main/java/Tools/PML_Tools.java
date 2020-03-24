package Tools;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.GaussianBlur3D;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Object3D_IJUtils;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Objects3DPopulationColocalisation;
import mcib3d.geom.PairColocalisation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import mcib3d.spatial.analysis.SpatialStatistics;
import mcib3d.spatial.descriptors.F_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;
        
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class PML_Tools {
   
    // Max distance to track object
    public static double maxDist = 5;
    // Background intensity factor for filter pml 
    public static double intFactor = 1.5; 
    public static boolean watershed = false; 
    
 /**
 * 
 * @param FileResults
 * @param resultsFileName
 * @param header
 * @return 
 */
public static BufferedWriter writeHeaders(String outDirResults, String resultsFileName, String header) throws IOException {
    FileWriter FileResults = new FileWriter(outDirResults + resultsFileName, false);
    BufferedWriter outPutResults = new BufferedWriter(FileResults); 
    outPutResults.write(header);
    outPutResults.flush();
    return(outPutResults);
}    
    
    /**
     * Dialog 
     */
    public static String dialog() {
        String dir = "";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addDirectoryField("Choose Directory Containing Image Files : ", "");
        gd.addNumericField("Threshold above diffuse PML intensity : ", intFactor, 2);
        gd.addCheckbox(" WaterShed split", watershed);
        gd.showDialog();
        dir = gd.getNextString();
        intFactor = gd.getNextNumber();
        watershed = gd.getNextBoolean();
        return(dir);
    }
    
    /**
     * Dialog when Zstep not found in metaData
     * @return zDepth
     */
    public static double dialogUnits(double zStep) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addMessage("Unable to find Z step !!");
        gd.addNumericField("Z step in µm : ", zStep, 3);
        gd.showDialog();
        double zDepth = gd.getNextNumber();
        if (gd.wasCanceled())
            zDepth = 0;
        return(zDepth);
    }
    
    /**Gaussian filter 
     * 
     * @param img
     * @param size
     */ 
    public static void gs_filter(ImagePlus img, double size) {
        GaussianBlur gaussian = new GaussianBlur();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            gaussian.blurGaussian(img.getProcessor(), size, size, 0.02);
            img.updateAndDraw();
        }
    }
    
    
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    
    /**
     * Nucleus segmentation for czi
     * @param imgNuc
     * @param radius
     * @return 
     */
    public static Objects3DPopulation find_nucleusCZI(ImagePlus imgNuc, int radius) {
        Calibration cal = imgNuc.getCalibration();
        IJ.run(imgNuc, "Remove Outliers...", "radius=20 threshold=1 which=Bright stack");
        IJ.run(imgNuc, "Difference of Gaussians", "  sigma1=30 sigma2=20 stack");
        threshold(imgNuc, AutoThresholder.Method.Li, true, false);
        Objects3DPopulation nucPop = new Objects3DPopulation();
        if (watershed) {
            ImagePlus imgWater = WatershedSplit(imgNuc, radius);
            imgWater.setCalibration(cal);
            nucPop = new Objects3DPopulation(imgWater);
            flush_close(imgWater);
        }
        else {
            nucPop = getPopFromImage(imgNuc, cal);
        }
        nucPop.removeObjectsTouchingBorders(imgNuc, false);
        return(nucPop);
    }
    
     /**
     * Nucleus segmentation
     * @param imgNuc
     * @param blur1
     * @param blur2
     * @param method
     * @param minNuc
     * @param maxNuc
     * @param radius
     * @param waterRad
     * @return 
     */
    public static Objects3DPopulation find_nucleus(ImagePlus imgNuc, String method, int blur1, int blur2, 
            int radius, int minNuc, int maxNuc, int waterRad) {
        Calibration cal = imgNuc.getCalibration();
        ImageStack stack = new ImageStack(imgNuc.getWidth(), imgNuc.getHeight());
        for (int i = 1; i <= imgNuc.getStackSize(); i++) {
            imgNuc.setZ(i);
            imgNuc.updateAndDraw();
            IJ.run(imgNuc, "Nuclei Outline", "blur="+blur1+" blur2="+blur2+" threshold_method="+method+
                    " outlier_radius="+radius+" outlier_threshold=1 max_nucleus_size="+maxNuc+
                    " min_nucleus_size="+minNuc+" erosion=5 expansion_inner=5 expansion=5 results_overlay");
            imgNuc.setZ(1);
            imgNuc.updateAndDraw();
            ImagePlus mask = new ImagePlus("mask", imgNuc.createRoiMask().getBufferedImage());
            ImageProcessor ip =  mask.getProcessor();
            ip.invertLut();
            for (int n = 0; n < 3; n++) 
                ip.erode();
            stack.addSlice(ip);
        }
        imgNuc.deleteRoi();
        ImagePlus imgStack = new ImagePlus("Nucleus", stack);
        imgStack.setCalibration(cal);
        Objects3DPopulation nucPop = new Objects3DPopulation();
        if (watershed) {
            ImagePlus imgWater = WatershedSplit(imgStack, radius);
            flush_close(imgStack);
            imgWater.setCalibration(cal);
            nucPop = new Objects3DPopulation(imgWater);
            flush_close(imgWater);
        }
        else {
            imgStack.setCalibration(imgNuc.getCalibration());
            nucPop = getPopFromImage(imgStack, cal);
        }
        nucPop.removeObjectsTouchingBorders(imgStack, false);
        flush_close(imgStack);
        return(nucPop);
    }
    
    /**
     * Size filter objects
     * remove touching border
     * 
     * @param min
     * @param max
     * @param objects
     * @param img
     * @param border
    */
    public static void objectsSizeFilter(double min, double max, Objects3DPopulation objects, ImagePlus img, boolean border) {
        ImageHandler imh = ImageInt.wrap(img).createSameDimensions();
        Object3D obj = null;
        boolean remove = false;
        if (objects.getNbObjects() > 0) {
            for (int n = 0; n < objects.getNbObjects(); n++) {
                remove = false;
                obj = objects.getObject(n);
                double vol = obj.getVolumeUnit();
                // remove if touching border
                if (border) {
                    if (obj.touchBorders(imh, false)) {
                        remove = true;
                    }
                }
                // Size filter remove small
                if ((vol < min) || (vol > max)) {
                    remove = true;
                }
                if (remove) {
                    objects.removeObject(n);
                    n--;
                }
            }
        }
    }
    
    // Threshold images and fill holes
    public static void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill, boolean calcul) {
        //  Threshold and binarize
       String cal = "";
       img.setZ(img.getNSlices()/2);
       img.updateAndDraw();
       IJ.setAutoThreshold(img, thMed.toString()+" dark");
       Prefs.blackBackground = false;
       if (calcul)
           cal = " calculate";
        IJ.run(img, "Convert to Mask", "method="+thMed.toString()+cal+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    public static Objects3DPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    /**
     * Read diffus dots intensity
     * fill dots voxel with zero in dots channel
     * Add to nucObj comment mean intensity pre-processing(integrated = false)
     * or integrated intensity (integrated = true)
     * 
     * @param dotsPop
     * @param nucObj
     * @param imgDots
     * @param integrated
     */
    public static void dotsDiffuse(Objects3DPopulation dotsPop, Nucleus nuc, ImagePlus imgDots, boolean integrated) {
        ImageHandler imhDotsDiffuse = ImageHandler.wrap(imgDots.duplicate());
        double dotsIntDiffuse ;
        int dots = dotsPop.getNbObjects();
        double volDotsDilated = 0;
        float dilate = 1.5f;
        for (int p = 0; p < dots; p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            // dilate 
            Object3DVoxels dotsDilatedObj = dotsObj.getDilatedObject(dilate, dilate, dilate);
            dotsDilatedObj.draw(imhDotsDiffuse, 0);
            volDotsDilated += dotsDilatedObj.getVolumeUnit();
        }
        double nucVolume = nuc.getVol();
        if (integrated)
            dotsIntDiffuse = (nuc.getObj().getIntegratedDensity(imhDotsDiffuse)/nucVolume) * (nucVolume - volDotsDilated); 
        else {
            ArrayUtil pixelInt = nuc.getObj().listValues(imhDotsDiffuse, 0);
            dotsIntDiffuse = pixelInt.getMean();
        }
        // put in nucleus object diffus intensity 
        nuc.setDiffuse(dotsIntDiffuse);
        nuc.setTotalInt(nuc.getObj().getIntegratedDensity(ImageHandler.wrap(imgDots)));
        imhDotsDiffuse.closeImagePlus();

    }
    
    /**
     * Create diffuse dots image
     * Fill zero in dots dots
     * @param dotsPop
     * @param nucObj
     * @param imgDotsOrg
     * @param outDirResults
     * @param rootName
     * @param diffuse
     * @param index
     * @param seriesName
     */
    public static void saveDiffuseImage(Objects3DPopulation dotsPop, Object3D nucObj, ImagePlus imgDotsOrg, String outDirResults,
            String rootName, String seriesName, String diffuse, int index) {
        ImageHandler imhDotsDiffuse = ImageHandler.wrap(imgDotsOrg.duplicate());
        String nucNumber = String.format("%03d", index);
        float dilate = 1.5f;
        for (int p = 0; p < dotsPop.getNbObjects(); p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            // dilate 
            Object3DVoxels dotsDilatedObj = dotsObj.getDilatedObject(dilate, dilate, dilate);
            dotsDilatedObj.draw(imhDotsDiffuse, 0);
        }
        ImagePlus imgColor = imhDotsDiffuse.getImagePlus();
        IJ.run(imgColor, "RGB Color", "");
        drawCountours(nucObj, imgColor, Color.white);
        // Save diffus
        FileSaver imgDiffus = new FileSaver(imgColor);
        imgDiffus.saveAsTiff(outDirResults + rootName + "_" + seriesName + "-Nuc"+nucNumber+"_"+diffuse+".tif");
        flush_close(imgColor);
        imhDotsDiffuse.closeImagePlus();
    }
    
    /**
     * Intensity filter objects
     * 
     * @param nucObj
     * @param dotsPop
     * @param imgDotsOrg
    */
    public static void ObjectsIntFilter(Nucleus nuc, Objects3DPopulation dotsPop, ImagePlus imgDotsOrg) {
        ImageHandler imhDotsOrg = ImageHandler.wrap(imgDotsOrg.duplicate());
        double dotsDiffuse = nuc.getDiffuse();
        // Remove dots with intensity < dots diffuse
        for (int p = 0; p < dotsPop.getNbObjects(); p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            double dotsInt = dotsObj.getPixMeanValue(imhDotsOrg);
            if (dotsInt <= dotsDiffuse*intFactor) {
                dotsPop.removeObject(p);
                p--;
            }
        }
        flush_close(imhDotsOrg.getImagePlus());
    }
    
    /*
    * Find PML inside nucleus
    * tag nucleus name as nucleus index
    * tag nucleus value as number of dots inside
    */
    public static Objects3DPopulation coloc(Object3D nucObj, Objects3DPopulation dotsPop) {
        Objects3DPopulation dotsInNuc = new Objects3DPopulation();
        for (int p = 0; p < dotsPop.getNbObjects(); p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            if (dotsObj.pcColoc(nucObj) > 80) {
                dotsInNuc.addObject(dotsObj);
            }
        }
        return(dotsInNuc);
    }
    
    /**
     * Get % of volume of coloc dost1 and dots2
     * @param nuc
     * @param dots1
     * @param dots2 
     */
    public static void findColoc(Nucleus nuc, Objects3DPopulation dots1, Objects3DPopulation dots2) {
        Objects3DPopulationColocalisation coloc = new Objects3DPopulationColocalisation(dots1, dots2);
        ArrayList<PairColocalisation> pairColoc = coloc.getAllColocalisationPairs();
        int volColocDots1 = pairColoc.get(0).getVolumeColoc();
        int volColocDots2 = pairColoc.get(1).getVolumeColoc();
        nuc.setColocDots1(volColocDots1/100);
        nuc.setColocDots2(volColocDots2/100);
    }
    
    
    // write object labels
    public static void labelsObject (Object3D obj, ImagePlus img, int number, int color) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, 32);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number), x, y);
        img.updateAndDraw();    
    }
    
    // tag object number with random color
    public static void tagsObject(ImageHandler imh, Object3D obj) {        
        int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
        obj.draw(imh, col);  
    }
    
   /**
     * ramdom color nucleus population
     */
    public static ImagePlus randomColorPop (Objects3DPopulation cellsPop,  ImageHandler img, boolean label) {
        //create image objects population
        img.set332RGBLut();
        img.setCalibration(img.getCalibration());
        for (int i = 0; i < cellsPop.getNbObjects(); i++) {
            Object3D obj = cellsPop.getObject(i);
            int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
            obj.draw(img, col);
            if (label)
               labelsObject(obj, img.getImagePlus(), (i+1), col); 
        } 
        return(img.getImagePlus());
    } 

    
    public static Point3D[] createEvaluationPoints(int numPoints, Objects3DPopulation population, Object3D mask) {
        Point3D[] evaluationPoints = new Point3D[numPoints];
        population.setMask(mask);
        for (int i = 0; i < numPoints; ++i) {
            evaluationPoints[i] = population.getRandomPointInMask();
        }
        return evaluationPoints;
    }
   
// Flush and close images
    public static void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
    
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public static void clearOutSide(ImagePlus img, Roi roi) {
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.updateAndDraw();
    }
    
    public static Plot createPlot(ArrayUtil xEvals, ArrayUtil[] sampleDistances, ArrayUtil observedDistances, ArrayUtil observedCD, ArrayUtil averageCD, String function) {
     
        Color ColorAVG = Color.blue;
        Color ColorENV = Color.gray;
        Color ColorOBS = Color.red;
        double plotMaxX = observedDistances.getMaximum();
        double plotMaxY = observedCD.getMaximum();
        int nbBins = 100;
        double env = 0.5;
        // low env
        double max = xEvals.getMaximum();
        ArrayUtil xEval0 = new ArrayUtil(nbBins);
        for (int i = 0; i < nbBins; i++) {
            xEval0.addValue(i, ((double) i) * max / ((double) nbBins));
        }
        // get the values
        ArrayUtil samplesPc5 = CDFTools.cdfPercentage(sampleDistances, xEval0, env / 2.0);
        ArrayUtil samplesPc95 = CDFTools.cdfPercentage(sampleDistances, xEval0, 1.0 - env / 2.0);
        // get the limits
        if (xEval0.getMaximum() > plotMaxX) {
            plotMaxX = xEval0.getMaximum();
        }
        if (samplesPc5.getMaximum() > plotMaxY) {
            plotMaxY = samplesPc5.getMaximum();
        }
        if (samplesPc95.getMaximum() > plotMaxY) {
            plotMaxY = samplesPc95.getMaximum();
        }
        if (xEvals.getMaximum() > plotMaxX) {
            plotMaxX = xEvals.getMaximum();
        }
        if (averageCD.getMaximum() > plotMaxY) {
            plotMaxY = averageCD.getMaximum();
        }
        if (observedCD.getMaximum() > plotMaxY) {
            plotMaxY = observedCD.getMaximum();
        }
        if (observedDistances.getMaximum() > plotMaxX) {
            plotMaxX = observedDistances.getMaximum();
        }
        // create the plot
        Plot plot = new Plot(function + "-function", "distance", "cumulated frequency");
        plot.setLimits(0, plotMaxX, 0, plotMaxY);

        // envelope  for e.g 10 % at 5 and 95 %
        plot.setColor(ColorENV);
        plot.addPoints(xEval0.getArray(), samplesPc5.getArray(), Plot.LINE);

        // envelope  for e.g 10 % at 5 and 95 %
        plot.setColor(ColorENV);
        plot.addPoints(xEval0.getArray(), samplesPc95.getArray(), Plot.LINE);

        // average
        plot.setColor(ColorAVG);
        plot.addPoints(xEvals.getArray(), averageCD.getArray(), Plot.LINE);

        // observed
        plot.setColor(ColorOBS);
        plot.addPoints(observedDistances.getArray(), observedCD.getArray(), Plot.LINE);

        return plot;
    }
    
    
     
    /**
    * For each nucleus compute F function
     * @param pop
     * @param mask
     * @param nuc
     * @param imgName
     * @param outDirResults
     * @return F SDI
    **/ 
    public static double processF (Objects3DPopulation pop, Object3D mask, String imgName, int nuc, String outDirResults) {

        // define spatial descriptor, model
        SpatialDescriptor spatialDesc = new F_Function(pop.getNbObjects(), mask);
        pop.setMask(mask);
        SpatialModel spatialModel = new SpatialRandomHardCore(pop.getNbObjects(), 0.8, mask);
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, 50, pop);
        spatialStatistics.setEnvelope(0.25);
        spatialStatistics.setVerbose(false);
//        Plot fPlot = spatialStatistics.getPlot();
//        fPlot.draw();
//        fPlot.addLabel(0.1, 0.1, "p = " + String.valueOf(spatialStatistics.getSdi()));
//        ImagePlus imgPlot = fPlot.getImagePlus();
//        FileSaver plotSave = new FileSaver(imgPlot);
//        plotSave.saveAsTiff(outDirResults + imgName + "_Fplot_" + nuc + ".tif");
//        flush_close(imgPlot);
        System.out.println("Nucleus" + nuc + " Sdi = " + spatialStatistics.getSdi());
        return(spatialStatistics.getSdi());
    }
    
    
    public static ImagePlus WatershedSplit(ImagePlus binaryMask, float rad) {
        float resXY = 1;
        float resZ = 1;
        float radXY = rad;
        float radZ = rad;
        Calibration cal = binaryMask.getCalibration();
        if (cal != null) {
            resXY = (float) cal.pixelWidth;
            resZ = (float) cal.pixelDepth;
            radZ = radXY * (resXY / resZ);
        }
        ImageInt imgMask = ImageInt.wrap(binaryMask);
        ImageFloat edt = EDT.run(imgMask, 0, resXY, resZ, false, 0);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 2.0, 2.0, 2.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        Watershed3D water = new Watershed3D(edt16, seedsImg, 0, 0);
        water.setLabelSeeds(true);
        return(water.getWatershedImage3D().getImagePlus());
    }
    
    
    
    /**
    * Compute dots results
    * Read intensity in mutant channel
    * @param nucObj nucleus
    * @param dotsPop dots population
    * @param imgDots read dots intensity
    * @param imgMut image file
     * @param results buffer
    **/
    public static void computeNucParameters(Object3D nucObj, int nucIndex, Objects3DPopulation dotsPop, ImagePlus imgDots, ImagePlus imgMut,
            String imgName, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");
        ImageHandler imhPML = ImageHandler.wrap(imgDots);
        ImageHandler imhMut;
        double  nucIntMut = 0;
        imhMut = ImageHandler.wrap(imgMut);
        nucIntMut = nucObj.getIntegratedDensity(imhMut);
        // measure nucleus volume
        // measure dots integrated intensity and volume
        double  nucVolume = nucObj.getVolumeUnit();
        double minDistBorder = Double.NaN;
        for (int p = 0; p < dotsPop.getNbObjects(); p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            double dotsIntensity = dotsObj.getIntegratedDensity(imhPML);
            double dotsVolume = dotsObj.getVolumeUnit();
            //double minDistBorder = dotsObj.distCenterBorderUnit(nucObj);
            results.write(imgName+"\t"+nucIndex+"\t"+nucVolume+"\t"+nucIntMut+"\t"+(p+1)+"\t"+dotsIntensity+"\t"+dotsVolume+"\t"+minDistBorder+"\n");
            results.flush();
        }
    }
    
     /**
    * Compute dots results
    * @param nucObj nucleus
    * @param dotsPop dots population
    * @param imgDots read dots intensity
     * @param imgName
     * @param results buffer
     * @throws java.io.IOException
    **/
    public static void computeNucParameters(Nucleus nuc, Objects3DPopulation dotsPop, ImagePlus imgDots,
            String imgName, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");
        ImageHandler imhPML = ImageHandler.wrap(imgDots);
        int dots = dotsPop.getNbObjects();
        double minDistBorder = Double.NaN;
        // measure nucleus volume
        // measure dots integrated intensity and volume
        for (int p = 0; p < dots; p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            double dotsIntensity = dotsObj.getIntegratedDensity(imhPML);
            double dotsVolume = dotsObj.getVolumeUnit();
            //minDistBorder = dotsObj.distCenterBorderUnit(nucObj);
            results.write(imgName+"\t"+nuc.getIndex()+"\t"+nuc.getVol()+"\t"+(p+1)+"\t"+dotsIntensity+"\t"+dotsVolume+"\t"+minDistBorder+"\n");
            results.flush();
        }
    }
    
    /**
    * Compute global nucleus and pml parameters for fixed cells
    * @param nucObj nucleus
     * @param nucIndex
    * @param pmlPop pml population
    * @param imgPML read pml intensity
     * @param imgMut
    * @param imgName image file
    * @param outDirResults results file
     * @param results buffer
     * @throws java.io.IOException
    **/
    public static void computeNucParameters2(Object3D nucObj, int nucIndex, Objects3DPopulation pmlPop, ImagePlus imgPML, ImagePlus imgMut,
            String imgName, String outDirResults, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");
        ImageHandler imhPML = ImageHandler.wrap(imgPML);
        ImageHandler imhMut;
        double  nucIntMut = 0;
        if (imgMut != null) {
            imhMut = ImageHandler.wrap(imgMut);
            nucIntMut = nucObj.getIntegratedDensity(imhMut);
        }
        // measure nucleus volume and integrated intensity in PML diffuse and image Mut
        // measure pml integrated intensity and volume
        DescriptiveStatistics pmlIntensity = new DescriptiveStatistics();
        DescriptiveStatistics pmlVolume = new DescriptiveStatistics();
        DescriptiveStatistics minDistBorder = new DescriptiveStatistics();
        double minDistCenterMean = Double.NaN;
        double minDistCenterSD = Double.NaN;
        double sdiF = Double.NaN;
        int pmlNuc = pmlPop.getNbObjects();
        double nucVolume = nucObj.getVolumeUnit();
        String nucIntDiffuse = nucObj.getComment();
        double nucShericity = nucObj.getSphericity(true);
        for (int p = 0; p < pmlNuc; p++) {
            Object3D pmlObj = pmlPop.getObject(p);
            pmlIntensity.addValue(pmlObj.getIntegratedDensity(imhPML));
            pmlVolume.addValue(pmlObj.getVolumeUnit());
            minDistBorder.addValue(pmlObj.distCenterBorderUnit(nucObj));
        }
//        if (pmlPop.getNbObjects() > 4) {
//            sdiF = processF(pmlPop, nucObj, imgName, nucIndex, outDirResults);
//        }
        if (pmlPop.getNbObjects() > 2) {
            minDistCenterMean = pmlPop.distancesAllClosestCenter().getMean(); 
            minDistCenterSD = pmlPop.distancesAllClosestCenter().getStdDev();
        }
        // compute statistics
        double pmlIntMean = pmlIntensity.getMean();
        double pmlIntSD = pmlIntensity.getStandardDeviation();
        double pmlIntMin = pmlIntensity.getMin();
        double pmlIntMax = pmlIntensity.getMax();
        double pmlVolumeMean = pmlVolume.getMean();
        double pmlVolumeSD = pmlVolume.getStandardDeviation();
        double pmlVolumeMin = pmlVolume.getMin();
        double pmlVolumeMax = pmlVolume.getMax();
        double pmlVolumeSum = pmlVolume.getSum();
        double minDistBorderMean = minDistBorder.getMean();
        double minDistBorderSD = minDistBorder.getStandardDeviation();

        results.write(imgName+"\t"+nucIndex+"\t"+nucVolume+"\t"+nucIntMut+"\t"+nucShericity+"\t"+pmlNuc+"\t"+nucIntDiffuse+"\t"+pmlIntMean+"\t"+
                pmlIntSD+"\t"+pmlIntMin+"\t"+pmlIntMax+"\t"+pmlVolumeMean+"\t"+pmlVolumeSD+"\t"+pmlVolumeMin+"\t"+pmlVolumeMax+"\t"+pmlVolumeSum+"\t"+
                minDistCenterMean+"\t"+minDistCenterSD+"\t"+minDistBorderMean+"\t"+minDistBorderSD+"\n");
        results.flush();
    }
    
    
    
     /**
    * Compute individual nucleus/dots parameters for live cells at time t
    * @param nucObj nucleus object
    * @param pmlPop pml population
    * @param imgPML read pml intensity
    * @param imgPMLDif read PML diffuse intensity
     * @param time
    * @param results buffer
    **/
    public static void computeNucParameters(Object3D nucObj, Objects3DPopulation pmlPop, ImagePlus imgPML, ImageHandler imgPMLDif, int time, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");

        // find dots inside nucleus
        // measure nucleus volume and integrated intensity in PML diffuse
        // measure pml integrated intensity and volume
            DescriptiveStatistics pmlIntensity = new DescriptiveStatistics();
            DescriptiveStatistics pmlVolume = new DescriptiveStatistics();
            int pmlinNuc = pmlPop.getNbObjects();
            ImageHandler imh = ImageHandler.wrap(imgPML);
            double nucVolume = nucObj.getVolumeUnit();
            double nucIntDiffuse = nucObj.getIntegratedDensity(imgPMLDif);
            double minDistCenterMean = pmlPop.distancesAllClosestCenter().getMean(); 
            double minDistCenterSD = pmlPop.distancesAllClosestCenter().getStdDev();
            int pmlIndex= 0;
            for (int p = 0; p < pmlPop.getNbObjects(); p++) {
                pmlIndex++;
                Object3D pml = pmlPop.getObject(p);
                pml.setName(String.valueOf(pmlIndex));
                pmlIntensity.addValue(pml.getIntegratedDensity(imh));
                pmlVolume.addValue(pml.getVolumeUnit());
            }
            // compute statistics
            double pmlIntMean = pmlIntensity.getMean();
            double pmlIntSD = pmlIntensity.getStandardDeviation();
            double pmlIntMin = pmlIntensity.getMin();
            double pmlIntMax= pmlIntensity.getMax();
            double pmlVolumeMean = pmlVolume.getMean();
            double pmlVolumeSD = pmlVolume.getStandardDeviation();
            double pmlVolumeMin = pmlVolume.getMin();
            double pmlVolumeMax = pmlVolume.getMax();
            double pmlVolumeSum = pmlVolume.getSum();
            
            results.write(time+"\t"+nucVolume+"\t"+pmlinNuc+"\t"+nucIntDiffuse+"\t"+pmlIntMean+"\t"+pmlIntSD+"\t"+pmlIntMin+"\t"+pmlIntMax+"\t"+
                    pmlVolumeMean+"\t"+pmlVolumeSD+"\t"+pmlVolumeMin+"\t"+pmlVolumeMax+"\t"+pmlVolumeSum+ "\t"+minDistCenterMean+"\t"+minDistCenterSD+"\n");
            results.flush();
    }
    
    
    /**
    * Compute global nucleus dots parameters for fixed cells
     * @param nuc
     * @param dots
    * @param dotsPop dots population
    * @param imgDots read dots intensity
     * @param imgName
     * @param results buffer
     * @throws java.io.IOException
    **/
    public static void computeGlobalNucParameters(Nucleus nuc, String dots, Objects3DPopulation dotsPop, ImagePlus imgDots, 
            String imgName, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");
        ImageHandler imhPML = ImageHandler.wrap(imgDots);
        // measure nucleus volume and integrated intensity in dots image
        // measure dots integrated intensity and volume
        DescriptiveStatistics dotsIntensity = new DescriptiveStatistics();
        DescriptiveStatistics dotsVolume = new DescriptiveStatistics();
        double dotsMinDistCenterMean = Double.NaN;
        double dotsMinDistCenterSD = Double.NaN;
        int dotsNuc = 0;
        double dotsVolColoc = 0;
        if (dots.equals("pml")) {
                dotsNuc = nuc.getDots1();
                dotsVolColoc = nuc.getColocDots1();
        }
        else {
             dotsNuc = nuc.getDots2();
             dotsVolColoc = nuc.getColocDots2();
        }
        double nucVolume = nuc.getVol();
        double nucShericity = nuc.getSphericity();
        double nucTotalPMLInt = nuc.getTotalInt();
        
        for (int p = 0; p < dotsNuc; p++) {
            Object3D dotsObj = dotsPop.getObject(p);
            dotsIntensity.addValue(dotsObj.getIntegratedDensity(imhPML));
            dotsVolume.addValue(dotsObj.getVolumeUnit());
        }

        if (dotsPop.getNbObjects() > 2) {
            dotsMinDistCenterMean = dotsPop.distancesAllClosestCenter().getMean(); 
            dotsMinDistCenterSD = dotsPop.distancesAllClosestCenter().getStdDev();
        }
        // compute statistics
        double dotsIntMean = dotsIntensity.getMean();
        double dotsIntSD = dotsIntensity.getStandardDeviation();
        double dotsIntMin = dotsIntensity.getMin();
        double dotsIntMax = dotsIntensity.getMax();
        double dotsVolumeMean = dotsVolume.getMean();
        double dotsVolumeSD = dotsVolume.getStandardDeviation();
        double dotsVolumeMin = dotsVolume.getMin();
        double dotsVolumeMax = dotsVolume.getMax();
        double dotsVolumeSum = dotsVolume.getSum();
        results.write(imgName+"\t"+nuc.getIndex()+"\t"+nucVolume+"\t"+nucShericity+"\t"+dotsNuc+"\t"+nuc.getDiffuse()+"\t"+
                nucTotalPMLInt+"\t"+dotsIntMean+"\t"+dotsIntSD+"\t"+dotsIntMin+"\t"+dotsIntMax+"\t"+dotsVolumeMean+"\t"+
                dotsVolumeSD+"\t"+dotsVolumeMin+"\t"+dotsVolumeMax+"\t"+dotsVolumeSum+"\t"+dotsMinDistCenterMean+"\t"+
                dotsMinDistCenterSD+"\t"+dotsVolColoc+"\n");
        results.flush();
    }
    
    
    /* create xml file for fixed cells
    * Write to xml file nucleus and pml parameters
    *   @param xmlFile file to write
    *   @param nucPop nucleus Population
    *   @param pmlPop pml population
    *   @param ims plm Diffuse image
    */

    public static void writeXml(String xmlFile, Objects3DPopulation nucPop, Objects3DPopulation pmlPop, ImageHandler imgPml, ImageHandler imgDiffus) throws ParserConfigurationException, TransformerConfigurationException, TransformerException {
        // create new XML doc
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        Document doc = docBuilder.newDocument();
        // XML root = nucleus
        Element rootElement = doc.createElement("Nucleus");
        doc.appendChild(rootElement);
        for (int n = 0; n < nucPop.getNbObjects(); n++) {
            Object3D nucObj = nucPop.getObject(n);
            // Add nucleus parameters
            Element nuc_param = doc.createElement("nucleus");
            rootElement.appendChild(nuc_param);
            // id
            Attr nucleus_id = doc.createAttribute("id");
            nucleus_id.setValue(nucObj.getName());
            nuc_param.setAttributeNode(nucleus_id);
            // volume
            Attr nucleus_vol = doc.createAttribute("volume");
            nucleus_vol.setValue(Double.toString(nucObj.getVolumeUnit()));
            nuc_param.setAttributeNode(nucleus_vol);
            // pml number
            Attr pml_nb = doc.createAttribute("pml_number");
            pml_nb.setValue(Integer.toString(nucObj.getValue()));
            nuc_param.setAttributeNode(pml_nb);
            // nucleus integrated intensity in diffuse pml image
            Attr diff_Int = doc.createAttribute("Diffuse_intensity");
            diff_Int.setValue(Double.toString(nucObj.getIntegratedDensity(imgDiffus)));
            nuc_param.setAttributeNode(diff_Int);
            for (int p = 0; p < pmlPop.getNbObjects(); p++) {
                Object3D pmlObj =  pmlPop.getObject(p);
                // write pml parameters for nucleus index
                if(pmlObj.getValue() == Integer.parseInt(nucObj.getName())) {
                    // Add pml parameters
                    Element pml_param = doc.createElement("pml");
                    nuc_param.appendChild(pml_param);
                    // id
                    Attr pml_id = doc.createAttribute("id");
                    pml_id.setValue(pmlObj.getComment());
                    pml_param.setAttributeNode(pml_id);
                    // Volume
                    Attr pml_vol = doc.createAttribute("volume");
                    pml_vol.setValue(Double.toString(pmlObj.getVolumeUnit()));
                    pml_param.setAttributeNode(pml_vol);
                    // Intensity
                    Attr pml_intensity = doc.createAttribute("intensity");
                    pml_intensity.setValue(Double.toString(pmlObj.getIntegratedDensity(imgPml)));
                    pml_param.setAttributeNode(pml_intensity);
                    // pml centroids           
                    Attr pml_XCentroid = doc.createAttribute("x_centroid");
                    pml_XCentroid.setValue(Double.toString(pmlObj.getCenterX()));
                    pml_param.setAttributeNode(pml_XCentroid);
                    Attr pml_YCentroid = doc.createAttribute("y_centroid");
                    pml_YCentroid.setValue(Double.toString(pmlObj.getCenterY()));
                    pml_param.setAttributeNode(pml_YCentroid);
                    Attr pml_ZCentroid = doc.createAttribute("z_centroid");
                    pml_ZCentroid.setValue(Double.toString(pmlObj.getCenterZ()));
                    pml_param.setAttributeNode(pml_ZCentroid);
                }
            }
        }
        // write the content into xml file
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(xmlFile));
        transformer.transform(source, result);
    }
    
    /* create xml file for live cells
    * Write to xml file nucleus and pml parameters
    *   @param xmlFile file to write
    *   @param nucPop nucleus Population
    *   @param pmlPop pml population
    *   @param ims plm Diffuse image
    */

    public static void writeXml(String xmlFile, Object3D nucObj, Objects3DPopulation pmlPop, ImageHandler imgPml, ImageHandler imgDiffus, int t) throws ParserConfigurationException, TransformerConfigurationException,
            TransformerException, SAXException, IOException {
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        Document doc;
        Element rootElement = null;
        // if t = 1 create new XML doc
        // and root as nucleus
        if (t == 1) {
            doc = docBuilder.newDocument();
            // XML root = nucleus
            rootElement = doc.createElement("Nucleus");
            doc.appendChild(rootElement);
            
        }
        // else parse xml
        else {
            doc = docBuilder.parse(xmlFile);
            rootElement = doc.getDocumentElement();
        }
        // Add time
        Element time = doc.createElement("time");
        rootElement.appendChild(time);
        // Add time index
        Attr t_index = doc.createAttribute("index");
        t_index.setValue(Integer.toString(t));
        time.setAttributeNode(t_index);
        // Add nucleus
        Element nucleus = doc.createElement("nucleus");
        time.appendChild(nucleus);
        // nucleus id
        Attr nucleus_id = doc.createAttribute("id");
        nucleus_id.setValue(nucObj.getName());
        nucleus.setAttributeNode(nucleus_id);
        // volume
        Attr nucleus_vol = doc.createAttribute("volume");
        nucleus_vol.setValue(Double.toString(nucObj.getVolumeUnit()));
        nucleus.setAttributeNode(nucleus_vol);
        // pml number
        Attr pml_nb = doc.createAttribute("pml_number");
        pml_nb.setValue(Integer.toString(nucObj.getValue()));
        nucleus.setAttributeNode(pml_nb);
        // nucleus integrated intensity in diffuse pml image
        Attr intensity = doc.createAttribute("diffuse_intensity");
        intensity.setValue(Double.toString(nucObj.getIntegratedDensity(imgDiffus)));
        nucleus.setAttributeNode(intensity);
        // add plm
        int pmlIndex = 0;
        for (int p = 0; p < pmlPop.getNbObjects(); p++) {
            pmlIndex++;
            Object3D pmlObj =  pmlPop.getObject(p);
            // Add pml parameters
            Element pml_param = doc.createElement("pml");
            nucleus.appendChild(pml_param);
            // id
            Attr pml_id = doc.createAttribute("id");
            pml_id.setValue(Integer.toString(pmlIndex));
            pml_param.setAttributeNode(pml_id);
            // Volume
            Attr pml_vol = doc.createAttribute("volume");
            pml_vol.setValue(Double.toString(pmlObj.getVolumeUnit()));
            pml_param.setAttributeNode(pml_vol);
            // Intensity
            Attr pml_intensity = doc.createAttribute("intensity");
            pml_intensity.setValue(Double.toString(pmlObj.getIntegratedDensity(imgPml)));
            pml_param.setAttributeNode(pml_intensity);
            // pml centroids           
            Attr pml_XCentroid = doc.createAttribute("x_centroid");
            pml_XCentroid.setValue(Double.toString(pmlObj.getCenterX()));
            pml_param.setAttributeNode(pml_XCentroid);
            Attr pml_YCentroid = doc.createAttribute("y_centroid");
            pml_YCentroid.setValue(Double.toString(pmlObj.getCenterY()));
            pml_param.setAttributeNode(pml_YCentroid);
            Attr pml_ZCentroid = doc.createAttribute("z_centroid");
            pml_ZCentroid.setValue(Double.toString(pmlObj.getCenterZ()));
            pml_param.setAttributeNode(pml_ZCentroid);
        }

        // write the content into xml file
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();
        DOMSource source = new DOMSource(doc);
        StreamResult result = new StreamResult(new File(xmlFile));
        transformer.transform(source, result);
    }

      /*
    Draw countours of objects
    */
    public static void drawCountours(Object3D obj, ImagePlus img, Color col) {
        ImagePlus imgMask = IJ.createImage("mask", img.getWidth(), img.getHeight(), img.getNSlices(), 8);
        for (int z = obj.getZmin(); z < obj.getZmax(); z++) {
            imgMask.setZ(z+1);
            ImageProcessor ip = imgMask.getProcessor();
            ByteProcessor bp = new ByteProcessor(ip, true);
            Object3D_IJUtils.draw(obj, bp, z, 255);
            ImagePlus maskPlus = new ImagePlus("mask " + z, bp);
            maskPlus.getProcessor().setAutoThreshold(AutoThresholder.Method.Default, true);
            IJ.run(maskPlus, "Create Selection", "");
            Roi roi = maskPlus.getRoi();
            img.setZ(z+1);
            img.getProcessor().setColor(col);
            img.getProcessor().drawRoi(roi);
            img.updateAndDraw();
        }   
        flush_close(imgMask);
    }
    
    
}
