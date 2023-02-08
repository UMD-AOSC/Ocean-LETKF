import yaml, sys

class Config:
    def __init__(self):
        pass

    def __str__(self):
        info="\n"
        for att in self.__dict__:
            info+="{} = {} \n".format(att,self.__dict__[att])
        return info

    def read_cfg(self, fnin, namelist="merge"):
        with open(fnin,"r") as f:
            cfg = yaml.safe_load(f)
        
        for att in cfg[namelist].keys(): 
            if hasattr(self, att):
               setattr(self, att, cfg[namelist][att])
            else:
               raise RuntimeError("Class ({}) does not have attribute ({})".format(type(self).__name__, att)) 
               sys.exit(146)
        print(self)
 
