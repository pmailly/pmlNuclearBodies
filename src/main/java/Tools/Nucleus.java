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
    private double totalInt;
    private int dots1;
    private int dots2;
    private double colocDots1;  // % of dots1 volume colocalize with dots2 population
    private double colocDots2;  // % of dots2 volume colocalize with dots1 population
   
	/**
         * 
         * @param nucObj
         * @param index
         * @param vol
         * @param sphericity
         * @param totalInt
         * @param diffuseInt
         * @param dots1 pml
         * @param dots2 dna or sumo 
         * @param colocDots1
         * @param colocDost2
         */
	public Nucleus(Object3D nucObj, int index, double vol, double sphericity, double totalInt, double diffuseInt,
                int dots1, int dots2, double colocDots1, double colocDots2) {
            this.nucObj = nucObj;
            this.index = index;
            this.vol = vol;
            this.sphericity = sphericity;
            this.diffuseInt = diffuseInt;
            this.dots1 = dots1;
            this.dots2 = dots2;
            this.colocDots1 = colocDots1;
            this.colocDots2 = colocDots2;
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
        
        public void setDiffuse(double diffuseInt) {
		this.diffuseInt = diffuseInt;
	}
        
        public void setTotalInt(double totalInt) {
		this.totalInt = totalInt;
	}
        
        public void setDots1(int dots1) {
		this.dots1 = dots1;
	}
        
        public void setDots2(int dots2) {
		this.dots2 = dots2;
	}
        
        public void setColocDots1(double colocDots1) {
		this.colocDots1 = colocDots1;
	}
        
        public void setColocDots2(double colocDots2) {
		this.colocDots2 = colocDots2;
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
            return(diffuseInt);
        }
        
        public double getTotalInt() {
            return(totalInt);
        }
        
        public int getDots1() {
            return dots1;
        }
        
	public int getDots2() {
		return dots2;
	}
        
        public double getColocDots1() {
		return colocDots1;
	}
        
        public double getColocDots2() {
		return colocDots2;
	}
}
