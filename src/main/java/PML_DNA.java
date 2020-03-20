/*
 * Find dots (PML) in nucleus DNA damages (Pierre)
 * version 2 for each nucleus object crop image before 
 * Measure integrated intensity, nb of dots per nucleus 
 * Author Philippe Mailly
 */


import static Tools.PML_Tools.writeHeaders;
import static Tools.PML_Tools.ObjectsIntFilter;
import static Tools.PML_Tools.coloc;
import static Tools.PML_Tools.computeGlobalNucParameters;
import static Tools.PML_Tools.computeNucParameters;
import static Tools.PML_Tools.dialog;
import static Tools.PML_Tools.dialogUnits;
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
import static Tools.PML_Tools.intFactor;
import static Tools.PML_Tools.labelsObject;
import static Tools.PML_Tools.median_filter;
import static Tools.PML_Tools.objectsSizeFilter;
import static Tools.PML_Tools.dotsDiffuse;
import static Tools.PML_Tools.randomColorPop;
import static Tools.PML_Tools.saveDiffuseImage;
import static Tools.PML_Tools.threshold;
import static Tools.PML_Tools.watershed;
import ij.gui.Roi;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import java.util.Arrays;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import org.apache.commons.io.FilenameUtils;


public class PML_DNA implements PlugIn {

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
            imageDir = dialog(false);
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
            String resultsName = "GlobalNucleusPMLResults_Int-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+".xls";
            String header = "ImageName\t#Nucleus\tNucleus Volume\tNucleus Sphericity\tPML dot number\tPML Diffuse IntDensity"
                    + "\tPML Mean dot IntDensity\tPML dot SD IntDensity\tPML dot Min IntDensity\tPML dot Max IntDensity"
                    + "\tPML dot Mean Volume\tPML dot SD Volume\tPML Min Vol\tPML Max Vol\tPML Sum Vol"
                    + "\tPML dot Mean center-center distance\tPML dot SD center-center distance"
                    + "\tPML dot Mean center-Nucleus border distance\tPML dot SD center-Nucleus border distance\n";
            BufferedWriter outPutPMLResultsGlobal = writeHeaders(outDirResults, resultsName, header); 
            
            // Global file for DNA results
            resultsName = "GlobalNucleusDNAResults_Int-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\tNucleus Sphericity\tDNA dot number\tDNA Diffuse IntDensity"
                    + "\tDNA Mean dot IntDensity\tDNA dot SD IntDensity\tDNA dot Min IntDensity\tDNA dot Max IntDensity"
                    + "\tDNA dot Mean Volume\tDNA dot SD Volume\tDNA Min Vol\tDNA Max Vol\tDNA Sum Vol"
                    + "\tDNA dot Mean center-center distance\tDNA dot SD center-center distance"
                    + "\tDNA dot Mean center-Nucleus border distance\tDNA dot SD center-Nucleus border distance\n";
            BufferedWriter outPutDNAResultsGlobal = writeHeaders(outDirResults, resultsName, header); 

            // Detailled parameters for PML results
            resultsName = "DetailledNucleusPMLResults_Int-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\t#PML dot\tPML dot IntDensity\tPML dot Volume"
                    + "\tPML dot center-center distance\n";
            BufferedWriter outPutPMLResultsDetail = writeHeaders(outDirResults, resultsName, header);
            
            // Detailled parameters for DNA results
            resultsName = "DetailledNucleusDNAResults_Int-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\t#DNA dot\tDNA dot IntDensity\tDNA dot Volume\t"
                    + "DNA dot center-center distance\n";
            BufferedWriter outPutDNAResultsDetail = writeHeaders(outDirResults, resultsName, header);
            
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
            for (int i = 0; i < imageFile.length; i++) {
                String fileExt = FilenameUtils.getExtension(imageFile[i]);
                if (fileExt.equals("nd")) {
                    imageName = inDir+ File.separator+imageFile[i];
                    rootName = imageFile[i].replace(".nd", "");
                    imageNum++;
                    reader.setId(imageName);
                    // Check calibration
                    if (imageNum == 1) {
                        series = 0;
                        cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                        cal.pixelHeight = cal.pixelWidth;
                        if (meta.getPixelsPhysicalSizeZ(series) != null)
                            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                        else
                            cal.pixelDepth = dialogUnits(zStep);
                        if (cal.pixelDepth == 0)
                            return;
                        cal.setUnit("microns");
                        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
                    }

                    series = reader.getSeriesCount();  
                    for (int s = 0; s < series; s++) {
                        reader.setSeries(s);
                        //seriesName = meta.getImageName(s);

                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(true);
                        options.setQuiet(true);
                        options.setSeriesOn(s, true);

                        /*
                        * Open DAPI channel
                        */
                        int dapiCh = 0;                        
                        options.setCBegin(s, dapiCh);
                        options.setCEnd(s, dapiCh);
                        System.out.println("-- Series : "+ seriesName);
                        System.out.println("Opening Nucleus channel");
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
                        ImgNucFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "_Nucleus-Objects.tif");
                        imhNuc.closeImagePlus();
                        flush_close(imgNuc);

                        // Open Orginal PML channel to read dot intensity
                        int pmlCh = 1;
                        System.out.println("Opening PML original channel");
                        options.setCBegin(s, pmlCh);
                        options.setCEnd(s, pmlCh);
                        ImagePlus imgDotsOrg = BF.openImagePlus(options)[0];
                        
                        // For all nucleus crop image
                        // Find PML in nucleus and DNA fragmentation
                        int nucIndex = 0;
                        
                        // Open DNA channel 
                        int dnaCh = 2;
                        System.out.println("Opening dna channel");
                        options.setCBegin(s, dnaCh);
                        options.setCEnd(s, dnaCh);
                        ImagePlus imgDNAOrg = BF.openImagePlus(options)[0];
                        
                        for (int n = 0; n < totalNucPop; n++) {
                            Object3D nucObj = nucPop.getObject(n);
                            int ZStartNuc = nucObj.getZmin() +1;
                            int ZStopNuc = nucObj.getZmax() + 1;
                            nucIndex++;
                            int[] box = nucObj.getBoundingBox();
                            Roi roiBox = new Roi(box[0], box[2], box[1] - box[0], box[3] - box[2]);
                            // Crop PML image
                            imgDotsOrg.setRoi(roiBox);
                            imgDotsOrg.updateAndDraw();
                            ImagePlus imgDotsCrop = new Duplicator().run(imgDotsOrg, ZStartNuc, ZStopNuc);
                            imgDotsCrop.deleteRoi();
                            imgDotsCrop.updateAndDraw();
                            ImagePlus imgDotsCropDup = imgDotsCrop.duplicate();
                            nucObj.translate(-nucObj.getXmin(), -nucObj.getYmin(), -ZStartNuc + 1);
                            nucObj.setName(Integer.toString(nucIndex));
                            
                            // Crop DNA image
                            imgDNAOrg.setRoi(roiBox);
                            imgDNAOrg.updateAndDraw();
                            ImagePlus imgDNACrop = new Duplicator().run(imgDNAOrg, ZStartNuc, ZStopNuc);
                            imgDNACrop.deleteRoi();
                            imgDNACrop.updateAndDraw();
                            ImagePlus imgDNACropDup = imgDNACrop.duplicate();
                            /*
                                Find pml dots in nucleus and compute diffuse intensity
                            */
                            median_filter(imgDotsCropDup, 1);
                            IJ.run(imgDotsCropDup, "Difference of Gaussians", " sigma1=4 sigma2=1 stack");
                            threshold(imgDotsCropDup, AutoThresholder.Method.RenyiEntropy, false, true);

                            Objects3DPopulation pmlPop = getPopFromImage(imgDotsCropDup, cal);
                            objectsSizeFilter(minPML, maxPML, pmlPop,imgDotsCropDup, false); 
                            System.out.println("PML pop after size filter = "+ pmlPop.getNbObjects());
                            // Find pml in nucleus
                            Objects3DPopulation pmlNucPop = coloc(nucObj, pmlPop);
                            System.out.println("Nucleus "+nucIndex+" PML = "+pmlNucPop.getNbObjects());
                            
                            // pre-processing PML diffus image intensity 
                            dotsDiffuse(pmlNucPop, nucObj, imgDotsCrop, false);
                             // intensity filter
                            //ObjectsIntFilter(nucObj, pmlNucPop, imgDotsCrop);
                            System.out.println("Nucleus "+nucIndex+" PML after intensity filter = "+pmlNucPop.getNbObjects());
                            // Find PML diffus intensity on pml filtered intensity
                            dotsDiffuse(pmlNucPop, nucObj, imgDotsCrop, true);
                            // save diffuse image
                            saveDiffuseImage(pmlNucPop, nucObj, imgDotsCrop, outDirResults, rootName, seriesName, "pmlDiffuse", nucIndex);
                            // Compute global PML parameters                        
                            // nucleus volume, nb of PML, mean PML intensity, mean PLM volume 
                            IJ.showStatus("Writing parameters ...");
                            computeGlobalNucParameters(nucObj, nucIndex, pmlNucPop, imgDotsCrop, rootName+seriesName+"_PML", outPutPMLResultsGlobal);
                            // Compute detalled PML parameters 
                            computeNucParameters(nucObj, pmlNucPop, imgDotsCrop, rootName+seriesName, outPutPMLResultsDetail);
                        
                            /*
                                Find dna fragmentation in nucleus and compute diffuse intensity
                            */
                            
                            IJ.run(imgDNACropDup, "Difference of Gaussians", " sigma1=4 sigma2=1 stack");
                            threshold(imgDNACropDup, AutoThresholder.Method.IJ_IsoData, false, true);

                            Objects3DPopulation dnaPop = getPopFromImage(imgDNACropDup, cal);
                            objectsSizeFilter(minPML, maxPML, dnaPop, imgDNACropDup, false); 
                            System.out.println("DNA pop after size filter = "+ dnaPop.getNbObjects());
                            // Find dna in nucleus
                            Objects3DPopulation dnaNucPop = coloc(nucObj, dnaPop);
                            System.out.println("Nucleus "+nucIndex+" DNA = "+dnaNucPop.getNbObjects());
                            
                            // pre-processing DNA diffus image intensity 
                            //dotsDiffuse(dnaNucPop, nucObj, imgDNACrop, false);
                             // intensity filter
                            //ObjectsIntFilter(nucObj, dnaNucPop, imgDNACrop);
                            //System.out.println("Nucleus "+nucIndex+" DNA after intensity filter = "+dnaNucPop.getNbObjects());
                            // Find DNA diffus intensity on DNA filtered intensity
                            dotsDiffuse(dnaNucPop, nucObj, imgDNACrop, true);
                            // save diffuse image
                            saveDiffuseImage(dnaNucPop, nucObj, imgDNACrop, outDirResults, rootName, seriesName, "DNA_Diffuse", nucIndex) ;
                            
                            // Compute global DNA parameters                        
                            // nucleus volume, nb of DNA, mean DNAL intensity, mean DNA volume 
                            IJ.showStatus("Writing parameters ...");
                            computeGlobalNucParameters(nucObj, nucIndex, dnaNucPop, imgDNACrop, rootName+seriesName, outPutDNAResultsGlobal);
                            // Compute detailled DNA parameters
                            computeNucParameters(nucObj,dnaNucPop, imgDNACrop, rootName+seriesName+"_DNA", outPutDNAResultsDetail);
                            
                            // Save objects image
                            ImageHandler imhPMLObjects = ImageHandler.wrap(imgDotsCrop).createSameDimensions();
                            ImageHandler imhNucObjects = imhPMLObjects.duplicate();
                            ImageHandler imhDNAObjects = imhPMLObjects.duplicate();
                            pmlNucPop.draw(imhPMLObjects, 255);
                            dnaNucPop.draw(imhDNAObjects, 255);
                            nucObj.draw(imhNucObjects, 255);
                            labelsObject(nucObj, imhNucObjects.getImagePlus(), nucIndex, 255);
                            ImagePlus[] imgColors = {imhDNAObjects.getImagePlus(), imhPMLObjects.getImagePlus(), imhNucObjects.getImagePlus()};
                            ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                            imgObjects.setCalibration(cal);
                            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                            FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                            ImgObjectsFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "-Nuc" + nucIndex + "-PML_Objects.tif");
                            flush_close(imgObjects);
                            flush_close(imhPMLObjects.getImagePlus());
                            flush_close(imhDNAObjects.getImagePlus());
                            flush_close(imhNucObjects.getImagePlus());
                            flush_close(imgDotsCropDup);
                            flush_close(imgDotsCrop);
                            flush_close(imgDNACrop);
                        }
                        flush_close(imgDotsOrg);
                        flush_close(imgDNAOrg);
                        options.setSeriesOn(s, false);
                    }
                }
            }
            try {
                outPutPMLResultsGlobal.close();
                outPutPMLResultsDetail.close();
                outPutDNAResultsGlobal.close();
                outPutDNAResultsDetail.close();
            } catch (IOException ex) {
                Logger.getLogger(PML_ES_NPM1C.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showStatus("Process done");
       
        }   catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(PML_ES_NPM1C.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}