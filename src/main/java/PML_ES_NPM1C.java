/*
 * Find dots (PML) in nucleus ES model NPM1C (Domitille)
 * version 2 for each nucleus object crop image before 
 * Measure integrated intensity, nb of dots per nucleus 
 * Author Philippe Mailly
 */

import static Tools.PML_Tools.segMethod;
import static Tools.PML_Tools.writeHeaders;
import static Tools.PML_Tools.ObjectsIntFilter;
import static Tools.PML_Tools.coloc;
import static Tools.PML_Tools.computeNucParameters;
import static Tools.PML_Tools.computeNucParameters2;
import static Tools.PML_Tools.dialog;
import static Tools.PML_Tools.dialogUnits;
import static Tools.PML_Tools.find_nucleusCZI;
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


public class PML_ES_NPM1C implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    public static final Calibration cal = new Calibration();

// min max volume in microns^3 for nucleus
     private final double minNuc = 50;
     private final double maxNuc = 5000;

// min volume in microns^3 for dots
    private final double minPML = 0.05;
// max volume in microns^3 for dots
    private final double maxPML = 10;

// Exist deconv image
    public boolean deconvImage = false;
// Default Z step
    private final double zStep = 0.153;

    

     

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
            imageDir = dialog(true);
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
             * Write headers results 
            */
            
            // Global parameters
            String resultsName = "GlobalNucleusResults_Int-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+".xls";
            String header = "ImageName\t#Nucleus\tNucleus Volume\tNucleus Int in Mut Channel\tNucleus Sphericity\tPML dot number"
                    + "\tPML Diffuse IntDensity\tPML Mean dot IntDensity\tPML dot SD IntDensity\tPML dot Min IntDensity"
                    + "\tPML dot Max IntDensity\tPML dot Mean Volume\tPML dot SD Volume\tPML Min Vol\tPML Max Vol"
                    + "\tPML Sum Vol\tPML dot Mean center-center distance\tPML dot SD center-center distance"
                    + "\tPML dot Mean center-Nucleus border distance\tPML dot SD center-Nucleus border distance\n";
            BufferedWriter outPutPMLResultsGlobal = writeHeaders(outDirResults, resultsName, header);
            
            // Detailled parameters
            resultsName = "DetailledPMLNucleusResults_Int-"+intFactor+"_WaterShed-"+Boolean.toString(watershed)+".xls";
            header = "ImageName\t#Nucleus\tNucleus Volume\tNucleus Int in Mut Channel\t#PML dot\tPML dot IntDensity"
                    + "\tPML dot Volume\tPML dot center-center distance\n";
            BufferedWriter outPutPMLResultsDetail = writeHeaders(outDirResults, resultsName, header);

            
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
                if (fileExt.equals("czi")) {
                    imageName = inDir+ File.separator+imageFile[i];
                    rootName = imageFile[i].replace(".czi", "");
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
                        if (zStep == 0)
                            return;
                        cal.setUnit("microns");
                        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
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
                        int dapiCh = 1;                        
                        options.setCBegin(s, dapiCh);
                        options.setCEnd(s, dapiCh);
                        System.out.println("-- Series : "+ seriesName);
                        System.out.println("Opening Nucleus channel");
                        ImagePlus imgNuc= BF.openImagePlus(options)[0];
                        Objects3DPopulation nucPop = new Objects3DPopulation();
                        if (segMethod == false)
                            nucPop = find_nucleus(imgNuc, "Li", 30, 50, 15, 20, 50, 5);
                        else
                            nucPop = find_nucleusCZI(imgNuc, 40);
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
                        int pmlCh = 0;
                        System.out.println("Opening PML original channel");
                        options.setCBegin(s, pmlCh);
                        options.setCEnd(s, pmlCh);
                        ImagePlus imgDotsOrg = BF.openImagePlus(options)[0];
                        
                        // For all nucleus crop image
                        // Find PML in nucleus
                        int nucIndex = 0;
                        // Open mut channel 
                        int mutCh = 2;
                        System.out.println("Opening red channel");
                        options.setCBegin(s, mutCh);
                        options.setCEnd(s, mutCh);
                        ImagePlus imgMutOrg = BF.openImagePlus(options)[0];
                        
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
                            // Crop mut image
                            imgMutOrg.setRoi(roiBox);
                            imgMutOrg.updateAndDraw();
                            ImagePlus imgMutCrop = new Duplicator().run(imgMutOrg, ZStartNuc, ZStopNuc);
                            imgMutCrop.deleteRoi();
                            imgMutCrop.updateAndDraw();
                            
                            // Detect pml dots
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
                            ObjectsIntFilter(nucObj, pmlNucPop, imgDotsCrop);
                            System.out.println("Nucleus "+nucIndex+" PML after intensity filter = "+pmlNucPop.getNbObjects());
                            // Find PML diffus intensity on pml filtered intensity
                            dotsDiffuse(pmlNucPop, nucObj, imgDotsCrop, true);
                            // save diffuse image
                            saveDiffuseImage(pmlNucPop, nucObj, imgDotsCrop, outDirResults, rootName, seriesName, "PML_Diffuse", nucIndex) ;
                            // Compute parameters                        
                            // nucleus volume, nb of PML, mean PML intensity, mean PLM volume
                            IJ.showStatus("Writing parameters ...");
                            computeNucParameters2(nucObj, nucIndex, pmlNucPop, imgDotsCrop, imgMutCrop, rootName+seriesName, outDirResults, outPutPMLResultsGlobal);
                            computeNucParameters(nucObj, pmlNucPop, imgDotsCrop, imgMutCrop, rootName+seriesName, outPutPMLResultsDetail);
                            
                            // Save objects image
                            ImageHandler imhDotsObjects = ImageHandler.wrap(imgDotsCrop).createSameDimensions();
                            ImageHandler imhNucObjects = imhDotsObjects.duplicate();
                            pmlNucPop.draw(imhDotsObjects, 255);
                            nucObj.draw(imhNucObjects, 255);
                            labelsObject(nucObj, imhNucObjects.getImagePlus(), nucIndex, 255);
                            ImagePlus[] imgColors = {imhDotsObjects.getImagePlus(), null, imhNucObjects.getImagePlus()};
                            ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                            imgObjects.setCalibration(cal);
                            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                            FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                            ImgObjectsFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "-Nuc" + nucIndex + "-PML_Objects.tif");
                            flush_close(imgObjects);
                            flush_close(imhDotsObjects.getImagePlus());
                            flush_close(imhNucObjects.getImagePlus());
                            flush_close(imgDotsCropDup);
                            flush_close(imgDotsCrop);
                            flush_close(imgMutCrop);
                        }
                        flush_close(imgDotsOrg);
                        flush_close(imgMutOrg);
                        options.setSeriesOn(s, false);
                    }
                }
            }
            try {
                outPutPMLResultsGlobal.close();
                 outPutPMLResultsDetail.close();
            } catch (IOException ex) {
                Logger.getLogger(PML_ES_NPM1C.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showStatus("Process done");
       
        }   catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(PML_ES_NPM1C.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}