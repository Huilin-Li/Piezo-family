import numpy as np
from . import tools


class TopologyCenters:
    def __init__(self, R, away, max_b, dff3_df, line3_df, max_a):
        self.start_center = (0, 0)
        self.away = int(away)
        self.max_a = int(max_a)
        self.max_b = int(max_b)
        self.R = float(R)
        self.df = dff3_df
        self.line3_df = line3_df
        self.centers = None
        self.TMCircleCenters_DICT = None
        self.membraney0 = None
        self.membraney1 = None
        self.Height = None
        self.TopoCenters = None
        self.TopoDataFrame = None

    def generate(self):
        df = self.df
        Nterm_IMO = df.iloc[0]["IMO"]
        if Nterm_IMO == "outside":
            self.OutsideNterm()
        else:
            self.InsideNterm()
        return self
    
    def AssignAminoAcide(self):
        TopoDataFrame = self.line3_df
        centers_list = self.TopoCenters
        print("centers_list", len(centers_list))
        TopoDataFrame["center"] = centers_list
        self.TopoDataFrame = TopoDataFrame
        return self
    


    def remove_duplicates(self, centers):
        CENTERS_arr = np.asarray(centers)
        CENTERS_arr = np.round(CENTERS_arr, 6)
        CENTERS_arr_list = list(map(tuple, CENTERS_arr.tolist()))
        seen = set()
        seen_add = seen.add
        return [x for x in CENTERS_arr_list if not (x in seen or seen_add(x))]
    
    def OutsideNterm(self):
        df = self.df
        for idx, row in df.iterrows():
            length = row["length"]
            if idx == len(df)-1:
                self.addCterm(length)
                break
            if idx == 0:
                # Nterm
                self.addNterm(length)
                self.membraney1 = self.centers[-1][-1]
                self.membraney0 = self.membraney1 - self.Height 
            elif idx%4 == 1:
                # download TM
                self.addDownwardingTMCenters(idx)
            elif idx%4 == 2:
                self.addIntracellularNotTMCenters(length) 
            elif idx%4 == 3:
                self.addUpwardingTMCenters(idx)
            else:
                self.addExtracellularNotTMCenters(length)

        self.TopoCenters = self.remove_duplicates(self.centers)
        self.AssignAminoAcide()
        return self
    
    def InsideNterm(self):
        df = self.df
        for idx, row in df.iterrows():
            length = row["length"]
            print("idx=", idx,"=", length)
            if idx == len(df)-1:
                self.addCterm(length)
                break
            if idx == 0:
                # Nterm
                self.addNterm(length)
                self.membraney0 = self.centers[-1][-1]
                self.membraney1 = self.membraney0 + self.Height 
            elif idx%4 == 1:
                # download TM
                self.addUpwardingTMCenters(idx) 
            elif idx%4 == 2:
                self.addExtracellularNotTMCenters(length) 
            elif idx%4 == 3:
                self.addDownwardingTMCenters(idx)
            else:
                self.addIntracellularNotTMCenters(length)

        self.TopoCenters = self.remove_duplicates(self.centers)
        self.AssignAminoAcide()
        return self

    def genTMCircleRelativeCenters(self, membraneThickness):
        TMUnits_idx_centers = tools.GenerateTMCircleRelativeCenters(df=self.df, membraneThickness=membraneThickness, R=self.R)
        self.TMCircleCenters_DICT = TMUnits_idx_centers
        TMUnits_idx_centers = self.TMCircleCenters_DICT
        Ybottoms = []
        Yups = []
        for k,v in TMUnits_idx_centers.items():
            Ybottoms.append(v[0][-1])
            Yups.append(v[-1][-1])
        bottom = sum(Ybottoms) / float(len(Ybottoms))
        up = sum(Yups) / float(len(Yups))
        self.Height = up - bottom
        return self
    

    def addNterm(self, length):
        df = self.df
        CENTERs_bridge_list = [self.start_center]
        Nterm_IMO = df.iloc[0]["IMO"]
        Nterm_centers_bridge = tools.AddNterm_Centers(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, IMO=Nterm_IMO, max_a=self.max_a, max_b=self.max_b)
        CENTERs_bridge_list += Nterm_centers_bridge[1:]
        print("add Nterm=CENTERs_bridge_list", len(CENTERs_bridge_list), len(set(CENTERs_bridge_list)))
        self.centers = CENTERs_bridge_list
        return self
    
    def addCterm(self, length):
        df = self.df
        CENTERs_bridge_list = self.centers
        Cterm_IMO = df.iloc[-1]["IMO"]
        Cterm_centers_bridge = tools.AddCterm_Centers(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, IMO=Cterm_IMO, max_a=self.max_a, max_b=self.max_b)
        CENTERs_bridge_list += Cterm_centers_bridge[1:]
        self.centers = CENTERs_bridge_list
        return self

    def addUpwardingTMCenters(self, idx):
        # upwarding
        CENTERs_bridge_list = self.centers
        upTMCenters_relative = self.TMCircleCenters_DICT[idx]
        df = self.df
        pre_IMO = df.iloc[idx-1]["IMO"]
        if pre_IMO == "outside":
            upTMCenters_relative = upTMCenters_relative[::-1]

        upTMCenters = tools.MoveCoords(coords=upTMCenters_relative, new_start=CENTERs_bridge_list[-1])
        CENTERs_bridge_list += upTMCenters
        print("idx=CENTERs_bridge_list", idx, "=", len(CENTERs_bridge_list), len(set(CENTERs_bridge_list)))
        self.centers = CENTERs_bridge_list
        return self
    
    def addDownwardingTMCenters(self, idx):
        # downwarding
        CENTERs_bridge_list = self.centers
        downTMCenters_relative = self.TMCircleCenters_DICT[idx]
        df = self.df
        pre_IMO = df.iloc[idx-1]["IMO"]
        if pre_IMO == "outside":
            downTMCenters_relative = downTMCenters_relative[::-1]
        
        downTMCenters = tools.MoveCoords(coords=downTMCenters_relative, new_start=CENTERs_bridge_list[-1])
        CENTERs_bridge_list += downTMCenters
        print("idx=CENTERs_bridge_list", idx, "=", len(CENTERs_bridge_list), len(set(CENTERs_bridge_list)))
        self.centers = CENTERs_bridge_list
        return self
    
    def addExtracellularNotTMCenters(self, length):
        CENTERs_bridge_list = self.centers
        ExtracellularNotTMCenters = tools.AddExtracellularNotTMCenters(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, max_a=self.max_a, max_b=self.max_b)
        CENTERs_bridge_list += ExtracellularNotTMCenters
        print("addExtracellularNotTMCenters=", len(CENTERs_bridge_list), len(set(CENTERs_bridge_list)))
        self.centers = CENTERs_bridge_list
        return self

    def addIntracellularNotTMCenters(self, length):
        CENTERs_bridge_list = self.centers
        IntracellularNotTMCenters = tools.AddIntracellularNotTMCenters(pre_center=CENTERs_bridge_list[-1], length=length, away=self.away, R=self.R, max_a=self.max_a, max_b=self.max_b)
        CENTERs_bridge_list += IntracellularNotTMCenters
        print("IntracellularNotTMCenters=", len(CENTERs_bridge_list), len(set(CENTERs_bridge_list)))
        self.centers = CENTERs_bridge_list
        return self
