# THIS PYTHON FILE CONTAINS AN ALGORITHM OF THE ANALYSIS AND DESIGN
# OF PILE FOUNDATIONS
# -------------------------------------------------------------------
class PileAnalysis:
    def __init__(
            self,
            support_df,
            pile_configuration,
            diameter,
            comp_cap,
            tens_cap,
            lat_cap,
            eta= 0.6,
            kappa= 0.4):
        self.support_df = support_df
        self.pile_configuration = pile_configuration
        self.diameter = diameter
        self.comp_cap = comp_cap
        self.tens_cap = tens_cap
        self.lat_cap = lat_cap
        self.eta = eta
        self.kappa = kappa

        self.__has_been_run = False

        # Numpy Library
        import numpy
        self.np = numpy

        # Pandas Library
        import pandas
        self.pd = pandas

    def run(self):
        from math import sqrt
        from statistics import mean
        # A. EFFICIENCY FACTORS
        # We will use a conservative method in computing the
        # efficiency factors for axial and lateral forces.
        
        # A1. Average distance between piles:
        distances = list()
        for p in self.pile_configuration:
            pile_conf = [
                pc for pc in self.pile_configuration
                if pc != p
            ]
            distances.append(
                min(
                    [sqrt( sum(
                        [ (p[0] - pc[0])**2, (p[1] - pc[1])**2 ]
                        ) )
                        for pc in pile_conf]
                    )
            )
        self.S = mean(distances)

        # A2. Axial Efficiency Factor:
        # The factor is computed in accordance with
        # (Sayed and Bakeer, 1992)
        eff_axial = 1 - (
            (1 - self.eta *self.kappa) *(10 / 7)
            *(self.tens_cap / self.comp_cap)
            )
        self.eff_axial = eff_axial

        # A3. Lateral Efficiency Factor
        # The factor is computed in accordance with
        # (Zhao and Stolarski, 1999)
        S_div_D = self.S / self.diameter
        eff_lat = 0.5 + 0.1 *S_div_D if 1 <= S_div_D < 5 \
                    else 1
        self.eff_lat = eff_lat

        # B. GENERAL ATTRIBUTE
        data_df = self.support_df.copy()
        pile_num = len(self.pile_configuration)
        data_len = len(data_df)
        Ix = sum( [pc[0]**2 for pc in self.pile_configuration] )
        Iy = sum( [pc[1]**2 for pc in self.pile_configuration] )
        x_max = max( [abs(pc[0]) for pc in self.pile_configuration] )
        y_max = max( [abs(pc[1]) for pc in self.pile_configuration] )
        
        # C. AXIAL FORCE (SINGLE PILE)
        single_comp = list()
        single_tens = list()
        Fz = list(data_df["Fz (kN)"])
        My = list(data_df["My (kNm)"])
        Mx = list(data_df["Mx (kNm)"])
        for a, b, c in zip(Fz, My, Mx):
            p_c = (a / pile_num) \
                + ( (abs(b) *x_max) / Ix ) \
                + ( (abs(c) *y_max) / Iy)
            single_comp.append(p_c)
            p_t = (a / pile_num) \
                - ( (abs(b) *x_max) / Ix ) \
                - ( (abs(c) *y_max) / Iy)
            if p_t < 0:
                single_tens.append(abs(p_t))
            else:
                single_tens.append(0)
        data_df["Pc (kN)"] = [float("%.2f" %r) for r in single_comp]
        data_df["Pt (kN)"] = [float("%.2f" %r) for r in single_tens]

        # D. LATERAL FORCE (SINGLE PILE)
        single_lat = list()
        Fx = list(data_df["Fx (kN)"])
        Fy = list(data_df["Fy (kN)"])
        for a, b in zip(Fx, Fy):
            h = sqrt(sum([a**2, b**2])) / pile_num
            single_lat.append(h)
        data_df["H (kN)"] = [float("%.2f" %r) for r in single_lat]

        # E. CAPACITY
        Pc_cap = [eff_axial *self.comp_cap] *data_len
        Pt_cap = [eff_axial *self.tens_cap] *data_len
        Hcap = [eff_lat *self.lat_cap] *data_len
        data_df["Pc_cap (kN)"] = [float("%.2f" %r) for r in Pc_cap]
        data_df["Pt_cap (kN)"] = [float("%.2f" %r) for r in Pt_cap]
        data_df["Hcap (kN)"] = [float("%.2f" %r) for r in Hcap]


        # F. RATIO
        # Ratio_comp = list()
        # Ratio_tens = list()
        # Ratio_lat = list()
        data_df["Ratio-C"] = data_df["Pc (kN)"] / data_df["Pc_cap (kN)"]
        data_df["Ratio-T"] = data_df["Pt (kN)"] / data_df["Pt_cap (kN)"]
        data_df["Ratio-L"] = data_df["H (kN)"] / data_df["Hcap (kN)"]

        data_df["Ratio-C"] = [float("%.3f" %r)
                                for r in data_df["Ratio-C"]]
        data_df["Ratio-T"] = [float("%.3f" %r)
                                for r in data_df["Ratio-T"]]
        data_df["Ratio-L"] = [float("%.3f" %r)
                                for r in data_df["Ratio-L"]]

        # G. SATISFIABILITY
        sat = list()
        for a, b, c in zip(
                data_df["Ratio-C"],
                data_df["Ratio-T"],
                data_df["Ratio-L"]
                ):
            sat.append("OK") if all(r < 1 for r in [a, b, c]) \
                else sat.append("NOT OK")
        # for k in range(data_len):
        #     if single_comp[k] <= Pc_cap[k] \
        #             and single_tens[k] <= Pt_cap[k] \
        #             and single_lat[k] <= Hcap[k]:
        #         sat.append("OK")
        #     else:
        #         sat.append("NOT OK")
        data_df["Sat."] = sat

        # H. INSTATIATING THE RESULT
        self.result = data_df.copy()

        # I. NOTIFICATION
        print(">>> The analysis has been performed []")
        self.__has_been_run = True

    def __not_run(self):
        print("ERROR: Please run the analysis first!")
        raise ValueError

    def present_result(self):
        if self.__has_been_run == True:
            display(self.result)
        else:
            self.__not_run()

    def save_excel(self, filename= "output.xlsx"):
        if self.__has_been_run == True:
            df = self.result.copy()
            df.to_excel(filename, index= False)
            print(">>> The result has been saved as an excel file []")
        else:
            self.__not_run()