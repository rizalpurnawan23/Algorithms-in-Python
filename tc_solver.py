class TowerCrane:
    def __init__(self,
                 pos_load_matrix,
                 pos_counter_matrix,
                 pos_weight_wind_tensor):
        self.pos_load_matrix = pos_load_matrix
        self.pos_counter_matrix = pos_counter_matrix
        self.pos_weight_wind_tensor = pos_weight_wind_tensor
        
        self.has_been_run = False
        
        # Numpy library:
        import numpy
        self.np = numpy
        
        # Pandas library:
        import pandas
        self.pd = pandas
                
    def read_me(self):
        print(
        """
        Class Name       : TowerCrane
        Type             : Python Class
        Date of Creation : 31/12/2022
        
        About the algorithm___
        
        This algorithm is created to compute the support reaction of
        a free standing tower crane. The underlying system used on
        the algorithm is modelled with linear algebra in accordance
        with classical static mechanics.
        
        The instantiations of the class are presented as follows:
        
        1. 'pos_load_matrix' consists of a list consisting of either
           two lists of 3 numbers or two numpy arrays of 3 numbers.
           The first entry in the list shall represent the position
           vector of the main crane's load. While the second entry
           shall represent the force vector of the main crane's load.
           
        2. 'pos_counter_vector' consists of a list consisting of
           either two lists of 3 numbers or two numpy arrays of 3
           numbers. The first entry in the list shall represent the
           position vector of the counterbalance load. While the
           second entry shall represent the force vector of the
           counter load.
           
        3. 'pos_weight_wind' consists of a list consisting of lists
           which consists of either lists of 3 numbers or numpy
           arrays of 3 numbers. For each list in the
           'pos_weight_wind', the first entry in the list shall
           represent a position vector in which a weight and
           a wind load apply, the second entry shall represent
           a force vector of the weight, and the third and the fourth
           entries shall represent a force vector of the wind load.
        
        --- Rizal Purnawan, the Author.
        """
        )
    
    def run(self):
        """
        Description___
        
        This function will run the algorithm such that we obtain
        the support reaction of the tower crane related to each
        applied force.
        """
        # A. Force Reactions:
        # A1. Due to Self Weight:
        sw_force_list = list()
        for w in self.pos_weight_wind_tensor:
            sw_force_list.append(self.np.array(w[1]))
        self.sw_force_reaction = -1 *sum(sw_force_list)
        
        # A2. Due to Main Loading:
        self.load_force_reaction = -1 *self.np.array(
            self.pos_load_matrix[1]
        )
        
        # A3. Due to Counterweight:
        self.counter_force_reaction = -1 *self.np.array(
            self.pos_counter_matrix[1]
        )
        
        # A4. Due to Wind Load:
        wind_force_list_X = list()      # X direction dominant
        wind_force_list_Y = list()      # Y direction dominant
        for w in self.pos_weight_wind_tensor:
            wind_force_list_X.append(self.np.array(w[2]))
            wind_force_list_Y.append(self.np.array(w[3]))
        self.wind_force_reaction_X = -1 *sum(wind_force_list_X)
        self.wind_force_reaction_Y = -1 *sum(wind_force_list_Y)
        
        # B. Moment Reactions
        # B1. Due to Self Weight:
        sw_moment_list = list()
        for w in self.pos_weight_wind_tensor:
            sw_moment_list.append(
                self.np.cross(w[0], w[1])
            )
        self.sw_moment_reaction = -1 *sum(sw_moment_list)
        
        # B2. Due to Main Loading:
        self.load_moment_reaction = -1 *self.np.cross(
            *self.pos_load_matrix
        )
        
        # B3. Due to Counterweight:
        self.counter_moment_reaction = -1 *self.np.cross(
            *self.pos_counter_matrix
        )
        
        # B4. Due to Wind Load:
        wind_moment_list_X = list()     # X direction dominant
        wind_moment_list_Y = list()     # Y direction dominant
        for w in self.pos_weight_wind_tensor:
            wind_moment_list_X.append(
                self.np.cross(w[0], w[2])
            )
            wind_moment_list_Y.append(
                self.np.cross(w[0], w[3])
            )
        self.wind_moment_reaction_X = -1 *sum(wind_moment_list_X)
        self.wind_moment_reaction_Y = -1 *sum(wind_moment_list_Y)
        
        # Setting the status of the class:
        self.has_been_run = True
        print(">>> The analysis has been completely performed []")

    def __not_run(self):
        print(">> ERROR: Please run the analysis first!")
        raise ValueError
        
    def force_reaction(self, comb= [1, 1, 1, 1], interactive= True):
        """
        Description___
        
        This function help us export the reaction force vector
        after the function 'run' has been executed. The parameter
        'comb' refers to the load combination being used. 'comb'
        shall be a list of 4 numbers. The detail description about
        'comb' is presented as follows:
        
        1. The first entry refers to the load combination for the
           self weight.
        2. The second entry refers to the load combination for the
           main loading of the tower crane.
        3. The third entry refers to the load combination for the
           counterweight.
        4. The fourth entry refers to the load combination for the
           wind load.

        On the otherhand, if the parameter 'interactive' is set to be
        'True', then this function returns an interactive output
        consisting of the information about the force reaction of the
        tower crane support. Otherwise, it returns a numpy array
        repersenting the vector of the support force reaction. By
        default, we set 'interactive' as 'True'.
        """
        if self.has_been_run == True:
            combined_force_list = [
                self.sw_force_reaction,
                self.load_force_reaction,
                self.counter_force_reaction,
                self.wind_force_reaction_X,
                self.wind_force_reaction_Y
            ]
            for k in range(len(combined_force_list)):
                combined_force_list[k] = comb[k] \
                                         *combined_force_list[k]
            force_reaction = sum(combined_force_list)
            if interactive == True:
                orientation = ["x", "y", "z"]
                print("Force Reaction:")
                for u, f in zip(orientation, force_reaction):
                    f = "%.3f" %f
                    f = float(f)
                    print(f"F{u} = {f} kN")
            else:
                return force_reaction
            
        else:
            self.__not_run()

    def moment_reaction(self, comb= [1, 1, 1, 1], interactive= True):
        """
        Description___

        This function help us export the reaction moment vector after
        the function 'run' has been executed. The parameter 'comb'
        refers to the load combination being used. 'comb' shall be a
        list of 4 numbers. The detail description of about 'comb' is
        presented as follows:

        1. The first entry refers to the load combination for the
           self weight.
        2. The second entry refers to the load combination for the
           main loading of the tower crane.
        3. The third entry refers to the load combination for the
           counterweight.
        4. The fourth entry refers to the load combination for the
           wind load.

        On the otherhand, if the parameter 'interactive' is set to be
        'True', then this function returns an interactive output
        consisting of the information about the moment reaction of the
        tower crane support. Otherwise, it returns a numpy array
        repersenting the vector of the support moment reaction. By
        default, we set 'interactive' as 'True'.
        """
        if self.has_been_run == True:
            combined_moment_list = [
                self.sw_moment_reaction,
                self.load_moment_reaction,
                self.counter_moment_reaction,
                self.wind_moment_reaction_X,
                self.wind_moment_reaction_Y
            ]
            for k in range(len(combined_moment_list)):
                combined_moment_list[k] = comb[k] \
                                          *combined_moment_list[k]
                moment_reaction = sum(combined_moment_list)
            if interactive == True:
                orientation = ["x", "y", "z"]
                print("Moment Reaction:")
                for u, m in zip(orientation, moment_reaction):
                    m = "%.3f" %m
                    m = float(m)
                    print(f"M{u} = {m} kNm")
            else:
                return moment_reaction
        else:
            self.__not_run()

    def comprehensive_output(self,
                             comb_list,
                             comb_name,
                             to_excel= [False, None],
                             show= True):
        """
        Description___

        This function is used to export a comprehensive output of the
        tower crane support.

        The parameter 'comb_list' shall be a list of load
        combinations in forms of lists of 4 numbers in which each of
        the numbers represents the load combination of the reaction
        element (self weight, main crane loading, counterweight,
        wind).

        The parameter 'to_excel' shall be a list of either a boolean
        value and a string or 'None'. If the boolean value is 'True',
        then the function exports the result into an excel data
        sheet. Otherwise it does not. The second entry is the excel
        file name.

        The parameter 'show' shall be either 'True' or 'False'. If
        it is set to be 'True', it displays the result in a form of
        pandas dataframe. Otherwise the function will return the
        pandas dataframe but not shown. By default we set the
        parameter to be 'True'.
        """
        if self.has_been_run == True:
            output_list = list()
            for lc, cn in zip(comb_list, comb_name):
                Fx, Fy, Fz = self.force_reaction(
                    comb= lc, interactive= False
                    )
                Mx, My, Mz = self.moment_reaction(
                    comb= lc, interactive= False
                    )
                all_forces = list()
                for f in [Fx, Fy, Fz, Mx, My, Mz]: # the order matters
                    f = "%.3f" %f
                    f = float(f)
                    all_forces.append(f)
                # output_list.append([cn, Fx, Fy, Fz, Mx, My, Mz])
                output_list.append([cn] + all_forces)
            columns = [
                "Comb",
                "Fx (kN)", "Fy (kN)", "Fz (kN)",
                "Mx (kNm)", "My (kNm)", "Mz (kNm)"
                ]
            output_df = self.pd.DataFrame(
                output_list, columns= columns
                )
            if to_excel[0] == True:
                file_name = to_excel[1] \
                            if not (to_excel[1] is None) \
                            else "tc_support_reaction"
                output_df.to_excel(file_name + ".xlsx", index= False)
            elif show == True:
                display(output_df)
            else:
                return output_df
        else:
            self.__not_run()