# Bismillahirrahmanirrahim
# Date of first creation: 05/28/2023

# DOCSTRING:
"""
    Title       : pile_foundation_solver
    Type        : Python Module
    Author      : Rizal Purnawan
    Date        : 05/28/2023

    About the module
    ----------------

    This module contains algorithms of either classes or functions
    or any executables for computing the capacity of piles in a
    group pile foundations as well as auxilliary attributes.

    Currently, this module contains the following python scripts:
    1. Class    -- PileFoundationSolver

    The readers should also read the technical documentation of each
    script which may contain the underlying mathematical description
    of the algorithm.

    -----------------------------------------------------------------
    Best regards,
    Rizal Purnawan, Author.
"""

# PILE FOUNDATION SOLVER:
class PileFoundationSolver:
    # ---------------------------------------------------------------
    # Docstring of the class:
    """
    Title       : PileFoundationSolver
    Type        : Python Class
    Author      : Rizal Purnawan
    Date        : 05/28/2023

    About the class
    ---------------

    This class is the main content of the script, which is a python
    class containing algorithms for computing the capacity of piles
    in a group pile foundation.

    The arguments for the class are described as follows:
    

    1. `pile_diameter`

        An integer or a floating point number representing the
        diameter of piles in use expressed in meters.


    2.  `column_coord`

        A list of two elements of numbers. Mathematically, this list
        represents the position vector of the structural column in
        the plane space of the pile cap. Each entry in the list
        represents the coordinate of the column in the plane.

        In general, the origin of the coordinates in the plane can be
        anywhere. However, for a convenient computational process,
        one can define the origin as the centre of the group pile.

        Note that each value in the list must take the unit of meter.

        
    3.  `list_of_piles_coord`
    
        A list containing lists of two elements of numbers. The
        elements of the largest list --which are lists, represent the
        position vectors of the piles in the plane space of the pile
        cap. The similar convention to the previous argument also
        applies for this argument.

        Note that each value in the lists must take the unit of
        meter.

        
    4.  `pile_cap_weight`

        A number (integer or floating point number) representing the
        weight of the pile cap expressed in kN.

        
    5.  `list_of_pile_capacity`

        Note that this algorithm does not compute the soil mechanical
        properties of the pile such as skin friction and end bearing.
        These soil mechanical properties must be prepared separately
        by the users before using this algorithm. And the required
        soil mechanical properties are the compressive capacity of
        the pile, the uplift/tensile capacity of the pile and the
        lateral capacity of the pile. All these values msut be the
        allowable value and must also be expressed in kN.

        The argument is a list containing the compressive, tensile
        and the lateral capacities of piles in the first, second and
        third entry of the list respectively.

        Note that all the numeric expression must be in kN.

    
    6.  `factors_of_safety`

        A list of two elements of numbers. The first entry must
        represent the factor of safety for the pile skin friction,
        and the second entry must represent the factor of safety
        for the pile end bearing.

    
    7. `sayed_bakeer_eta_kappa`

        A list of two elements of numbers. The first entry must
        represent the eta prime factor for the Sayed-Bakeer axial
        efficiency factor. While the second entry must represent the
        kappa factor for the same axial efficiency factor. By
        default, we define this argument as [0.6, 0.4].

        In the current developement, we only use axial efficiency
        factor of group pile proposed by (Sayed and Bakeer, 1992).
        On the other hand, we use the lateral efficiency factor
        proposed by (Zhao and Stolarski, 1999).

    
    8.  `additional_coords_loads`

        A list of two lists of three elements of numbers. The first
        entry in the largest list --which is a list of three
        numbers-- represents the position vector of an additional
        load or force on the pile cap. And the second entry of the
        largest list --which is a list of three numbers-- represents
        the force vector of the additional load or force on the pile
        cap. By default, we set it to be [[0, 0, 0], [0, 0, 0]].


    Additional libraries to be used in this class are Numpy and
    Pandas, which will be imported and attributed in the init
    function.
    """

    # ---------------------------------------------------------------
    # `init` function:
    def __init__(self,
                 pile_diameter,
                 column_coord,
                 list_of_piles_coord,
                 pile_cap_weight,
                 list_of_pile_capacity,
                 factors_of_safety,
                 sayed_bakeer_eta_kappa= [0.6, 0.4],
                 additional_coords_loads= [[0, 0, 0], [0, 0, 0]]):
        try:
            # All the necessary clauses are made here:
            if all(type(x) in [int, float]
                        for x in [pile_diameter, pile_cap_weight]
                        ) \
                    and all(type(li) == list and len(li) == 2
                            and all(type(x) in [int, float]
                                    for x in li)
                            for li in [column_coord,
                                       factors_of_safety,
                                       sayed_bakeer_eta_kappa]) \
                    and type(list_of_piles_coord) == list \
                    and all(type(li) == list and len(li) == 2
                            and all(type(x) in [int, float]
                                    for x in li)
                            for li in list_of_piles_coord) \
                    and type(additional_coords_loads) == list \
                    and all(type(li) == list and len(li) == 3
                            and all(type(x) in [int, float]
                                    for x in li)
                            for li in additional_coords_loads):
                pass
            else:
                raise ValueError
            
            # If the clauses are met, here are the instantations of
            # the attributes:
            self.pile_diameter = pile_diameter
            self.colum_coord = column_coord
            self.list_of_piles_coord = list_of_piles_coord
            self.pile_cap_weight = pile_cap_weight
            self.list_of_pile_capacity = list_of_pile_capacity
            self.factors_of_safety = factors_of_safety
            self.sayed_bakeer_eta_kappa = sayed_bakeer_eta_kappa
            self.additional_coords_loads = additional_coords_loads

            # Attributes for the output:
            self.result_df = None
            self.eff_axial = None
            self.eff_lateral = None

            # In addition, since we also use Numpy and Pandas, we
            # instatiate them into attributes as follows:
            import numpy as np
            import pandas as pd
            from math import sqrt
            from statistics import mean
            self.__np = np
            self.__pd = pd
            self.__sqrt = sqrt
            self.__mean = mean
        except:
            print("ERROR: Invalid arguments!")
            print(">>> Please read the docstring carefully and give valid arguments.")
            raise ValueError
        
    # ---------------------------------------------------------------
    # `analyse_pile` Function:
    def analyse_pile(self,
                     excel_path,
                     sheet_name= 0,
                     seismic= False,
                     tensile_to_skin_friction_ratio= 0.7,
                     output= "df"
                     ):
        """
        Title       : analyse_pile
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 05/28/2023

        About the function
        ------------------

        This function will execute the computation of the pile
        capacity in a group pile. The arguments given in the function
        are described as follows:

        1.  `excel_path` and `sheet_name`

            A path of an excel file containing the preliminary data
            for computing the capacity of a pile in a group.
            `excel_path` must be a string representing the path of
            the excel file, while `sheet_name` is a string
            representing the sheet name in the excel file. We have
            made a conventional format for this matter, and the users
            shall follow this convention in order to use the
            algorithm.
            
            The convention is described as follows:
            
            1a) In the current development, the units to be used are
                meters (m) for length and kilo Newton (kN) for forces
                and loads.
            
            1b) The data in the excel spread sheet of the preliminary
                data shall be of the following format:
                -----------------------------------------------------
                >   The first row shall be the title of the columns.
                -----------------------------------------------------
                >   Column 1 of the preliminary data shall contain
                    the load combinations being considered. Thus, the
                    data type of column 1 can either be string or
                    integer.
                -----------------------------------------------------
                >   Column 2 to column 4 shall be the forces of the
                    support reaction of the structural column from
                    the structural modelling. Column 2 and column 3
                    shall be the forces in the horizontal directions.
                    While column 4 shall be forces in the vertical
                    direction.
                    NOTE: The unit for the forces shall be kN.
                -----------------------------------------------------
                >   Column 5 to column 7 shal be the moments of the
                    support reaction of the structural column from
                    the structural modelling. Column 5 and column 6
                    shall be the moments with respect to the
                    horizontal axes of the spacial model. While
                    column 7 shall be the moments with respect to the
                    vertical axis.
                    NOTE: The unit for the moments shall be kNm.
                -----------------------------------------------------

        2.  `seismic`

            A boolean value of either True or False. By default it is
            set to be False. This argument indicate whether the
            special seismic consideration is taken into account. If
            it does, then the pile capacity will be multiplied by a
            factor of 1.2. Otherwise, the pile capcaity remains in
            its initial state.

        3.  `tensile_to_skin_friction_ratio`

            A number in between 0 and 1 representing the ratio
            between the tensile capacity of piles and the skin
            friction of piles. By default, we set the value to
            be 0.7, which means the tensile capcity is 70% of the
            skin friction in this consideration.

        4.  `output`

            A string of either 'complete' or 'df'. If 'complete',
            the function returns a triple of the dataframe result
            of the pile analysis, the axial efficiency factor and
            the lateral efficiency factor. If 'df', the function
            only retursn the dataframe result of the pile analysis.

        
        The function returns a data frame containing the result of
        the analysis which can be instantiated outside the class.
        """
        # -----------------------------------------------------------
        # Instantiating necessary properties/attributes:
        np = self.__np
        pd = self.__pd
        sqrt = self.__sqrt
        mean = self.__mean

        pile_diameter = self.pile_diameter
        column_coord = self.colum_coord
        list_of_piles_coord = self.list_of_piles_coord
        pile_cap_weight = self.pile_cap_weight
        # Here is how the seismic consideration takes place:
        if seismic == True:
            fac = 1.2
        else:
            fac = 1
        c_all = fac *self.list_of_pile_capacity[0]
        t_all = fac *self.list_of_pile_capacity[1]
        l_all = fac *self.list_of_pile_capacity[2]
        # Other instantiations:
        sf_f, sf_p = self.factors_of_safety
        eta_prime, kappa = self.sayed_bakeer_eta_kappa
        add_coord, add_force = self.additional_coords_loads

        # -----------------------------------------------------------
        # Converting the excel_path into a pandas dataframe:
        try:
            prelim_data = pd.read_excel(
                excel_path,
                sheet_name= sheet_name
                )
            df_columns = [
                "LC",
                "F1", "F2", "F3",
                "M1", "M2", "M3"
                ]
            prelim_data.columns = df_columns
            
        except:
            print("ERROR: Invalid directory!")
            print(">>> Please use a valid directory and a valid data.")
            raise ValueError
        
        # -----------------------------------------------------------
        # Computing the ultimate skin friction and end bearing:
        # Ultimate skin friction:
        Qf = sf_f *t_all / tensile_to_skin_friction_ratio
        # Ultimate end bearing:
        Qp = sf_p *(c_all - t_all / tensile_to_skin_friction_ratio)

        # -----------------------------------------------------------
        # Computing the axial efficiency factor:
        # As mentoined, in the current development, the axial safety
        # factor uses the one proposed by (Sayed and Bakeer, 1992):
        eff_axial = 1 - (1 - eta_prime *kappa) *(Qf / (Qf + Qp))

        # -----------------------------------------------------------
        # Computing the lateral efficiency factor:
        # Expected distance of adjacent piles:
        adj_dist = max(
            [
                min([np.linalg.norm(np.array(p) - np.array(x))
                for p in list_of_piles_coord if p != x])
                for x in list_of_piles_coord
                ]
        )
        if 1 < adj_dist / pile_diameter < 5:
            eff_lateral = 0.5  + 0.1 *(adj_dist / pile_diameter)
        elif adj_dist / pile_diameter >= 5:
            eff_lateral = 1

        # -----------------------------------------------------------
        # Computing the center of piles:
        pile_center = [
            # Computing the center of piles in direction 1 and
            # direction 2 iteratively:
            mean([p[i] for p in list_of_piles_coord])
            for i in [0, 1]
        ]

        # -----------------------------------------------------------
        # Computing the expected eccentricity of piles:
        # Direction 1:
        p1_max = max(
            [
                # Taking the absolute inner product of a unit vector
                # [1, 0] and each pile's position vector minus the
                # position vector of the center of group pile.
                abs(np.dot(
                        [1, 0], np.array(p) - np.array(pile_center))
                    )
                for p in list_of_piles_coord
            ]
        )
        p1_sum = sum(
            [
                # Taking the square of inner product of a unit vector
                # [1, 0] and each pile's position vector minus the
                # position vector of the center of group pile.
                np.dot(
                    [1, 0], np.array(p) - np.array(pile_center)
                )**2
                for p in list_of_piles_coord
                ]
        )
        # The following designation shall be made to avoid division
        # by zero:
        if p1_max == 0 or p1_sum == 0:
            p1_max, p1_sum = 0, 1
        # Direction 2:
        p2_max = max(
            [
                # Taking the absolute inner product of a unit vector
                # [0, 1] and each pile's position vector minus the
                # position vector of the center of group pile.
                abs(np.dot(
                        [0, 1], np.array(p) - np.array(pile_center))
                    )
                for p in list_of_piles_coord
            ]
        )
        p2_sum = sum(
            [
                # Taking the square of inner product of a unit vector
                # [0, 1] and each pile's position vector minus the
                # position vector of the center of group pile.
                np.dot(
                    [0, 1], np.array(p) - np.array(pile_center)
                )**2
                for p in list_of_piles_coord
                ]
        )
        # The following designation shall be made to avoid division
        # by zero:
        if p2_max == 0 or p2_sum == 0:
            p2_max, p2_sum = 0, 1

        # -----------------------------------------------------------
        # Executing iterations:
        list_of_output = list()
        for k in range(len(prelim_data)):
            # List of row:
            row = list()

            # Recapturing the load combinations:
            load_comb = prelim_data[df_columns[0]][k]
            row.append(load_comb)

            # Force Vector:
            # Note that column 2 to column 4 in the prelim_data
            # are the forces of the support reaction.
            vecF = np.array(
                list(prelim_data.loc[k, :])[1:4]
                )
            row = row + list(vecF)

            # Moment Vector:
            # Note that column 5 to column 7 in the prelim_data
            # are the moments of the support reaction.
            vecM = np.array(
                list(prelim_data.loc[k, :])[4:]
            )
            row = row + list(vecM)

            # Moment added with possible column eccentricity and
            # an additional load:
            vec_col = np.array(
                column_coord + [0]
            )
            pile_center_3D = np.array(
                pile_center + [0]
            )
            vecM = vecM \
                    + np.cross(vec_col - pile_center_3D, vecF) \
                    + np.cross(add_coord - pile_center_3D, add_force)
            
            # Computing maximum compression in a single pile:
            pile_compression = max(
                (
                    np.dot([0, 0, 1], vecF + [0, 0, pile_cap_weight])
                    / len(list_of_piles_coord)
                ) + (
                    abs(np.dot([0, 1, 0], vecM)) *p1_max / p1_sum
                ) + (
                    abs(np.dot([1, 0, 0], vecM)) *p2_max / p2_sum
                ),
                0
            )
            row.append(round(pile_compression, 3))

            # Computing maximum tension in a single pile:
            pile_tension = max(
                -(
                    (
                        np.dot([0, 0, 1],
                               vecF + [0, 0, pile_cap_weight]
                               )
                        / len(list_of_piles_coord)
                    ) - (
                        abs(np.dot([0, 1, 0], vecM)) *p1_max / p1_sum
                    ) - (
                        abs(np.dot([1, 0, 0], vecM)) *p2_max / p2_sum
                    )
                ),
                0
            )
            row.append(round(pile_tension, 3))

            # Computing maximum lateral force in a single pile:
            pile_lateral = abs(np.dot([1, 1, 0], vecF)) \
                            / (len(list_of_piles_coord) *sqrt(2))
            row.append(round(pile_lateral, 3))

            # Allowable compression:
            row.append(
                round(
                    eff_axial *c_all, 3
                )
            )

            # Allowable tension:
            row.append(
                round(
                    eff_axial *t_all, 3
                )
            )

            # Allowable lateral capacity:
            row.append(
                round(
                    eff_lateral *l_all, 3
                )
            )

            # Compression ratio:
            rc = pile_compression / (eff_axial *c_all)
            row.append(round(rc, 3))

            # Tension ratio:
            rt = pile_tension / (eff_axial *t_all)
            row.append(round(rt, 3))

            # Lateral capacity ratio:
            rl = pile_lateral / (eff_lateral *l_all)
            row.append(round(rl, 3))

            # Satisfiability notification:
            if all(r <= 1 for r in [rc, rt, rl]):
                row.append("SAT")
            else:
                row.append("UNSAT")

            # Final output row:
            list_of_output.append(row)

        # -----------------------------------------------------------
        # Forming a dataframe of the output:
        result_df = pd.DataFrame(list_of_output)
        result_df.columns = [
            "Load Comb",
            "F1 (kN)", "F2 (kN)", "F3 (kN)",
            "M1 (kNm)", "M2 (kNm)", "M3 (kNm)",
            "P_comp (kN)", "P_tens (kN)", "P_lat (kN)",
            "Eff_ax C_all (kN)",
            "Eff_ax T_all (kN)",
            "Eff_lat L_all",
            "rc", "rt", "rl",
            "SAT"
        ]

        # -----------------------------------------------------------
        # Instantiating results into class attributes: [CANCELLED]
        # self.result_df = result_df
        # self.eff_axial = eff_axial
        # self.eff_lateral = eff_lateral

        # -----------------------------------------------------------
        # Product of the function:
        if output == "complete":
            return result_df, eff_axial, eff_lateral
        elif output == "df":
            return result_df
        else:
            print("ERROR: Invalid 'output' argument!")
            raise ValueError
    
    # ---------------------------------------------------------------
    def save_result(self,
                    excel_path,
                    sheet_name= 0,
                    file_name= "output.xlsx",
                    column_names= None):
        """
        Title       : save_result
        Type        : Class Function
        Author      : Rizal Purnawan
        Date        : 05/29/2023

        About the function
        ------------------

        This function help save the result of the computation
        conducted by executing the class function `analyse_pile`.
        The arguments of the function are described comprehensively
        as follows:

        1.  `excel_path` and `sheet_name`

            Read the docstring for function `analyse_pile`.
        
        2.  `file_name`

            The excel directory of the saved file. By default it is
            set to be 'output.xlsx'.

        3.  `column_names`

            The users are free to modify the names of the columns of
            the resulting dataframe. By default it is set to be None,
            and it will use the column names from the initial
            definition.
        """
        # It must be made sure that the function `analyse_pile`
        # should have been run before users run this function.
        result_df, eff_axial, eff_lateral \
            = self.analyse_pile(excel_path,
                                sheet_name= sheet_name,
                                output= "complete")
        if column_names is not None:
            if type(column_names) == list \
                    and all(type(cn) == str for cn in column_names) \
                    and len(column_names) == 17:
                pass
            else:
                print("ERROR: Use a proper column names!")
                raise ValueError
        else:
            column_names = result_df.column

        # Saving the result:
        result_df.column = column_names
        result_df.to_excel(file_name, index= False)
        print(">>> Result saved []")