/*
 * Find  PML and SUMO dots in nucleus (Pierre)
 * version 2 for each nucleus object crop image before 
 * Measure integrated intensity, nb of dots per nucleus 
 * Author Philippe Mailly
 */

import Tools.Nucleus;
import static Tools.PML_Tools.writeHeaders;
import static Tools.PML_Tools.ObjectsIntFilter;
import static Tools.PML_Tools.coloc;
import static Tools.PML_Tools.computeGlobalNucParameters;
import static Tools.PML_Tools.computeNucParameters;
import static Tools.PML_Tools.find_nucleus;
import static Tools.PML_Tools.getPopFromImage;
import ij.*;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.geom.Object3D;
import static Tools.PML_Tools.flush_close;
import static Tools.PML_Tools.labelsObject;
import static Tools.PML_Tools.median_filter;
import static Tools.PML_Tools.objectsSizeFilter;
import static Tools.PML_Tools.dotsDiffuse;
import static Tools.PML_Tools.findColoc;
import static Tools.PML_Tools.intFactor;
import static Tools.PML_Tools.randomColorPop;
import static Tools.PML_Tools.saveDiffuseImage;
import static Tools.PML_Tools.threshold;
import static Tools.PML_Tools.watershed;
import fiji.util.gui.GenericDialogPlus;
import ij.gui.Roi;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;
import java.util.Arrays;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


public class PML_SUMO implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    public static final Calibration cal = new Calibration();

// min max volume in microns^3 for nucleus
     private final double minNuc = 20;
     private final double maxNuc = 2000;

// min volume in microns^3 for dots
    private final double minPML = 0.05;
// max volume in microns^3 for dots
    private final double maxPML = 10;

// Exist deconv image
    public boolean deconvImage = false;
   
// Default Z step
    public static double zStep = 0.193;
// Nucleus     
    public static Nucleus nucleus = new Nucleus(null, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    

    /**
     * Dialog ask for channels order
     * @param channels
     * @param showCal
     * @param cal
     * @return ch;

     */
    public ArrayList dialog(String[] channels, boolean showCal, Calibration cal) {
        ArrayList ch = new ArrayList();
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("DAPI       : ", channels, channels[0]);
        gd.addChoice("PML        : ", channels, channels[1]);
        gd.addChoice("SUMO       : ", channels, channels[2]);
        gd.addNumericField("Threshold above diffuse PML intensity : ", intFactor, 2);
        if (showCal) {
            gd.addMessage("No Z step calibration found", Font.getFont("Monospace"), Color.red);
            gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
            gd.addNumericField("Z pixel size : ", zStep, 3);
        }
        gd.showDialog();
        ch.add(0, gd.getNextChoice());
        ch.add(1, gd.getNextChoice());
        ch.add(2, gd.getNextChoice());
        if (showCal) {
            cal.pixelWidth = gd.getNextNumber();
            cal.pixelDepth = gd.getNextNumber();
        }
        
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing Image Files...");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Results_IntFactor-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            /** 
             * Write headers results for results file
             */
            // Global file for PML results
            String resultsName = "GlobalNucleusPMLResults_Int-"+intFactor+".xls";
            String header = "ImageName\t#Nucleus\tNucleus Volume\tNucleus Sphericity\tPML dot number\tPML Total IntDensity"
                    + "\tPML Diffuse IntDensity\tPML Mean dot IntDensity\tPML dot SD IntDensity\tPML dot Min IntDensity"
                    + "\tPML dot Max IntDensity\tPML dot Mean Volume\tPML dot SD Volume\tPML Min Vol\tPML Max Vol"
                    + "\tPML Sum Vol\tPML dot Mean center-center distance\tPML dot SD center-center distance\tPML Coloc\n";
            BufferedWriter outPutPMLResultsGlobal = writeHeaders(outDirResults, resultsName, header); 
            
            // Global file for SUMO results
            resultsName = "GlobalNucleusSUMOResults_Int-"+intFactor+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\tNucleus Sphericity\tSUMO dot number\tSUMO Total IntDensity"
                    + "\tSUMO Diffuse IntDensity\tSUMO Mean dot IntDensity\tSUMO dot SD IntDensity\tSUMO dot Min IntDensity"
                    + "\tSUMO dot Max IntDensity\tSUMO dot Mean Volume\tSUMO dot SD Volume\tSUMO Min Vol\tSUMO Max Vol"
                    + "\tSUMO Sum Vol\tSUMO dot Mean center-center distance\tSUMO dot SD center-center distance\tSUMO Coloc\n";
            BufferedWriter outPutSUMOResultsGlobal = writeHeaders(outDirResults, resultsName, header); 

            // Detailled parameters for PML results
            resultsName = "DetailledNucleusPMLResults_Int-"+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\t#PML dot\tPML dot IntDensity\tPML dot Volume"
                    + "\tPML dot center-center distance\n";
            BufferedWriter outPutPMLResultsDetail = writeHeaders(outDirResults, resultsName, header);
            
            // Detailled parameters for SUMO results
            resultsName = "DetailledNucleusSUMOResults_Int-"+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\t#SUMO dot\tSUMO dot IntDensity\tSUMO dot Volume\t"
                    + "SUMO dot center-center distance\n";
            BufferedWriter outPutSUMOResultsDetail = writeHeaders(outDirResults, resultsName, header);
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0;
            int series;
            String imageName = "";
            String rootName = "";
            String seriesName = "";
            ArrayList<String> ch = new ArrayList();
            for (int i = 0; i < imageFile.length; i++) {
                String fileExt = FilenameUtils.getExtension(imageFile[i]);
                if (fileExt.equals("czi")) {
                    imageName = inDir+ File.separator+imageFile[i];
                    rootName = imageFile[i].replace(".czi", "");
                    imageNum++;
                    boolean showCal = false;
                    reader.setId(imageName);
                    int chNb = reader.getSizeC();
                    String[] channels = new String[chNb];
                    for (int c = 0; c < chNb; c++) 
                        channels[c] = meta.getChannelFluor(0, c);
                    // Check calibration
                    if (imageNum == 1) {
                        series = 0;
                        cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                        cal.pixelHeight = cal.pixelWidth;
                        if (meta.getPixelsPhysicalSizeZ(series) != null)
                            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                        else
                            showCal= true;
                        if (zStep == 0)
                            return;
                        cal.setUnit("microns");
                        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
                        
                        // return the index for channels DAPI, Astro, Dots and ask for calibration if needed 
                        ch = dialog(channels, showCal, cal);

                        if (ch == null) {
                            IJ.showStatus("Plugin cancelled !!!");
                            return;
                        }
                    }

                    series = reader.getSeriesCount();  
                    for (int s = 0; s < series; s++) {
                        reader.setSeries(s);
                        seriesName = meta.getImageName(s);

                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(true);
                        options.setQuiet(true);
                        options.setSeriesOn(s, true);

                        /*
                        * Open DAPI channel
                        */
                        int channelIndex = ArrayUtils.indexOf(channels, ch.get(0));
                        options.setCBegin(s, channelIndex);
                        options.setCEnd(s, channelIndex);
                        System.out.println("-- Series : "+ seriesName);
                        System.out.println("Opening Nucleus channel "+ ch.get(0));
                        ImagePlus imgNuc= BF.openImagePlus(options)[0];
                        
                        Objects3DPopulation nucPop = new Objects3DPopulation();
                        nucPop = find_nucleus(imgNuc, "Triangle", 10, 30, 15, 5, 20, 8);
                        objectsSizeFilter(minNuc, maxNuc, nucPop, imgNuc, true);
                        int totalNucPop = nucPop.getNbObjects();
                        System.out.println("nucleus after size filter= "+totalNucPop);
                        
                        // tags nucleus
                        ImageHandler imhNuc = ImageInt.wrap(imgNuc).createSameDimensions();
                        randomColorPop(nucPop, imhNuc, true);
                        // save image for nucleus population
                        imhNuc.getImagePlus().setCalibration(cal);
                        FileSaver ImgNucFile = new FileSaver(imhNuc.getImagePlus());
                        ImgNucFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "_Nucleus_Objects.tif");
                        imhNuc.closeImagePlus();
                        flush_close(imgNuc);

                        // Open Original PML channel to read dot intensity
                        channelIndex = ArrayUtils.indexOf(channels, ch.get(1));
                        System.out.println("Opening PML original channel"+ch.get(1));
                        options.setCBegin(s, channelIndex);
                        options.setCEnd(s, channelIndex);
                        ImagePlus imgPMLOrg = BF.openImagePlus(options)[0];
                        
                        // For all nucleus crop image
                        // Find PML and SUMO dots in nucleus
                        int nucIndex = 0;
                        
                        // Open SUMO channel 
                        channelIndex = ArrayUtils.indexOf(channels, ch.get(2));
                        System.out.println("Opening SUMO channel"+ch.get(2));
                        options.setCBegin(s, channelIndex);
                        options.setCEnd(s, channelIndex);
                        ImagePlus imgSUMOOrg = BF.openImagePlus(options)[0];
                        
                        for (int n = 0; n < totalNucPop; n++) {
                            nucIndex++;
                            Object3D nucObj = nucPop.getObject(n);
                            
                            int ZStartNuc = nucObj.getZmin() +1;
                            int ZStopNuc = nucObj.getZmax() + 1;
                            
                            int[] box = nucObj.getBoundingBox();
                            Roi roiBox = new Roi(box[0], box[2], box[1] - box[0], box[3] - box[2]);
                            
                            // Crop PML image
                            imgPMLOrg.setRoi(roiBox);
                            imgPMLOrg.updateAndDraw();
                            ImagePlus imgPMLCrop = new Duplicator().run(imgPMLOrg, ZStartNuc, ZStopNuc);
                            imgPMLCrop.deleteRoi();
                            imgPMLCrop.updateAndDraw();
                            ImagePlus imgPMLCropDup = imgPMLCrop.duplicate();
                            nucObj.translate(-nucObj.getXmin(), -nucObj.getYmin(), -ZStartNuc + 1);
                            Nucleus nucleusObj = new Nucleus(nucObj, nucIndex, nucObj.getVolumeUnit(), nucObj.getSphericity(true),
                                    0, 0, 0, 0, 0, 0);
                            
                            /*
                                Find pml dots in nucleus and compute diffuse intensity
                            */
                            median_filter(imgPMLCropDup, 1);
                            IJ.run(imgPMLCropDup, "Difference of Gaussians", " sigma1=3 sigma2=1 stack");
                            threshold(imgPMLCropDup, AutoThresholder.Method.RenyiEntropy, false, false);

                            Objects3DPopulation pmlPop = getPopFromImage(imgPMLCropDup, cal);
                            objectsSizeFilter(minPML, maxPML, pmlPop,imgPMLCropDup, false); 
                            System.out.println("PML pop after size filter = "+ pmlPop.getNbObjects());
                            // Find pml in nucleus
                            Objects3DPopulation pmlNucPop = coloc(nucObj, pmlPop);
                            System.out.println("Nucleus "+nucIndex+" PML = "+pmlNucPop.getNbObjects());
                            
                            // pre-processing PML diffus image intensity 
                            dotsDiffuse(pmlNucPop, nucleusObj, imgPMLCrop, false);
                             // intensity filter
                            ObjectsIntFilter(nucleusObj, pmlNucPop, imgPMLCrop);
                            System.out.println("Nucleus "+nucIndex+" PML after intensity filter = "+pmlNucPop.getNbObjects());
                            
                            // Find PML diffus intensity on pml filtered intensity
                            dotsDiffuse(pmlNucPop, nucleusObj, imgPMLCrop, true);
                            // save diffuse image
                            //saveDiffuseImage(pmlNucPop, nucObj, imgPMLCrop, outDirResults, rootName, seriesName, "PML_Diffuse", nucIndex);
                            // add pml number to Nucleus
                            nucleusObj.setDots1(pmlNucPop.getNbObjects());
                            
                            // Compute detailed PML parameters 
                            computeNucParameters(nucleusObj, pmlNucPop, imgPMLCrop, rootName+seriesName, outPutPMLResultsDetail);
                        
                            /*
                                Find SUMO dots in nucleus and compute diffuse and total intensity
                            */
                            // Crop SUMO image
                            imgSUMOOrg.setRoi(roiBox);
                            imgSUMOOrg.updateAndDraw();
                            ImagePlus imgSUMOCrop = new Duplicator().run(imgSUMOOrg, ZStartNuc, ZStopNuc);
                            imgSUMOCrop.deleteRoi();
                            imgSUMOCrop.updateAndDraw();
                            ImagePlus imgSUMOCropDup = imgSUMOCrop.duplicate();
                            
                            IJ.run(imgSUMOCropDup, "Difference of Gaussians", " sigma1=3 sigma2=1 stack");
                            threshold(imgSUMOCropDup, AutoThresholder.Method.RenyiEntropy, false, false);

                            Objects3DPopulation SumoPop = getPopFromImage(imgSUMOCropDup, cal);
                            objectsSizeFilter(minPML, maxPML, SumoPop, imgSUMOCropDup, false); 
                            System.out.println("SUMO pop after size filter = "+ SumoPop.getNbObjects());
                            // Find Sumo in nucleus
                            Objects3DPopulation sumoNucPop = coloc(nucObj, SumoPop);
                            System.out.println("Nucleus "+nucIndex+" SUMO = "+sumoNucPop.getNbObjects());
                            
                            // pre-processing SUMO diffus image intensity 
                            dotsDiffuse(sumoNucPop, nucleusObj, imgSUMOCrop, false);
                             // intensity filter
                            ObjectsIntFilter(nucleusObj, sumoNucPop, imgSUMOCrop);
                            System.out.println("Nucleus "+nucIndex+" SUMO after intensity filter = "+sumoNucPop.getNbObjects());
                            // Find SUMO diffus intensity on SUMO filtered intensity
                            dotsDiffuse(sumoNucPop, nucleusObj, imgSUMOCrop, true);
                            // save diffuse image
                            //saveDiffuseImage(sumoNucPop, nucObj, imgSUMOCrop, outDirResults, rootName, seriesName, "SUMO_Diffuse", nucIndex) ;
                            // add Sumo number to Nucleus
                            nucleusObj.setDots2(sumoNucPop.getNbObjects());
                            
                            // Compute detailled DNA parameters
                            computeNucParameters(nucleusObj, sumoNucPop, imgSUMOCrop, rootName+seriesName+"_DNA", outPutSUMOResultsDetail);
                                          
                            // Find colocalization in PML and SUMO populations
                            findColoc(nucleusObj, pmlNucPop, sumoNucPop);
                            
                            // Compute global PML parameters                        
                            // nucleus volume, nb of PML, mean PML intensity, mean PLM volume 
                            IJ.showStatus("Writing parameters ...");
                            computeGlobalNucParameters(nucleusObj, "pml", pmlNucPop, imgPMLCrop, rootName+seriesName+"_PML", outPutPMLResultsGlobal);
                            
                            // Compute global DNA parameters                        
                            // nucleus volume, nb of DNA, mean DNAL intensity, mean DNA volume 
                            IJ.showStatus("Writing parameters ...");
                            computeGlobalNucParameters(nucleusObj, "sumo", sumoNucPop, imgSUMOCrop, rootName+seriesName, outPutSUMOResultsGlobal);
                            
                            // Save objects image
                            String nucNumber = String.format("%03d", nucIndex);
                            ImageHandler imhPMLObjects = ImageHandler.wrap(imgPMLCrop).createSameDimensions();
                            ImageHandler imhNucObjects = imhPMLObjects.duplicate();
                            ImageHandler imhSUMOObjects = imhPMLObjects.duplicate();
                            pmlNucPop.draw(imhPMLObjects, 255);
                            sumoNucPop.draw(imhSUMOObjects, 255);
                            nucObj.draw(imhNucObjects, 255);
                            labelsObject(nucObj, imhNucObjects.getImagePlus(), nucIndex, 255);
                            ImagePlus[] imgColors = {imhSUMOObjects.getImagePlus(), imhPMLObjects.getImagePlus(), imhNucObjects.getImagePlus()};
                            ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                            imgObjects.setCalibration(cal);
                            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                            FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                            ImgObjectsFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "-Nuc" + nucNumber + "_PML_Objects.tif");
                            flush_close(imgObjects);
                            flush_close(imhPMLObjects.getImagePlus());
                            flush_close(imhSUMOObjects.getImagePlus());
                            flush_close(imhNucObjects.getImagePlus());
                            flush_close(imgPMLCropDup);
                            flush_close(imgPMLCrop);
                            flush_close(imgSUMOCrop);
                        }
                        flush_close(imgPMLOrg);
                        flush_close(imgSUMOOrg);
                        options.setSeriesOn(s, false);
                    }
                }
            }
            try {
                outPutPMLResultsGlobal.close();
                outPutPMLResultsDetail.close();
                outPutSUMOResultsGlobal.close();
                outPutSUMOResultsDetail.close();
            } catch (IOException ex) {
                Logger.getLogger(PML_ES_NPM1C.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showStatus("Process done");
       
        }   catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(PML_ES_NPM1C.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}