/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tools;

import mcib3d.geom.Object3D;

/**
 *
 * @author phm
 */
public class Nucleus {
    private Object3D nucObj;
    private int index;
    private double vol;
    private double sphericity;
    private double diffuseInt;
    private int pml;
    private int dna;
   
	
	public Nucleus(Object3D nucObj, int index, double vol, double sphericity, double diffuseInt, int pml, int dna) {
            this.nucObj = nucObj;
            this.index = index;
            this.vol = vol;
            this.sphericity = sphericity;
            this.diffuseInt = diffuseInt;
            this.pml = pml;
            this.dna = dna;
	}
        
         public void setObj(Object3D nucObj) {
		this.nucObj = nucObj;
	}
         
        public void setIndex(int index) {
		this.index = index;
	}
        
        public void setVol(double vol) {
		this.vol = vol;
	}
        
        public void setSphericity(double sphericity) {
		this.sphericity = sphericity;
	}
        
        public void setDiffuseInt(double diffuseInt) {
		this.diffuseInt = diffuseInt;
	}
        
        public void setPML(int pml) {
		this.pml = pml;
	}
        
        public void setDNA(int dna) {
		this.dna = dna;
	}
        
        public Object3D getObj() {
            return nucObj;
        }
              
        public int getIndex() {
            return index;
        }
        
        public double getVol() {
            return(vol);
        }
        
        public double getSphericity() {
            return(sphericity);
        }
        
        public double getDiffuse() {
            return(sphericity);
        }
        
        public int getPML() {
            return pml;
        }
        
	public int getDNA() {
		return dna;
	}
        
}
