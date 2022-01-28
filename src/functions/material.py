# -*- coding: utf-8 -*-

class Material:
    def __init__(material_code, geometry):
        if (material_code == "2pu"):
            self.nu = np.array([(2.84)])
            self.nu = np.reshape(nu,(1,1))
            
            self.sigf = np.array([(0.0816)])
            self.sigf = np.reshape(sigf,(1,1))
            
            self.sigc = np.array([(0.019584)])
            self.sigc = np.reshape(sigc,(1,1))
            
            self.sigs = np.array([(0.225216)])
            self.sigs = np.reshape(sigs,(1,1,1))
            
            self.sigt = np.array([(0.32640)])
            self.sigt = np.reshape(sigt,(1,1))
            
            self.X = np.array([(1.0)])
            self.X = np.reshape(X,(1,1))
            
            #critical radius of geo
            if (Geo == 1):
                self.R = np.array([2.256751]) #slab
                self.true_flux = np.array([0.9701734, 0.8810540, 0.7318131, 0.4902592])
                self.true_flux_r =  np.array([0.25,0.5,0.75,1.0])*self.R[0]
            elif (Geo == 2):
                self.R = np.array([4.27996]) #cylinder
            else:
                self.R = np.array([6.082547]) #sphere
        else:
            print("Type 'help(Material)' for a list of available materials")
    
        