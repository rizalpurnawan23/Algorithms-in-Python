# Algoritma
class OptimasiLenturBalok:
    def readMe(self):
        return (
        """
        ___Tentang algoritma
        
            Algoritma ini merupakan suatu tool untuk melakukan analisis
        pada struktur lentur beton bertulang dengan mempertimbangkan
        optimasi hasil analisis. Hal yang sangat menguntungkan dalam
        penggunaan suatu bahasa program adalah kemampuan bahasa program
        untuk melakukan komputasi dengan lebih leluasa, terlebih dalam hal
        iterasi. Program ini diharapkan dapat menjadi manfaat untuk para
        praktisi dalam bidang structural engineering khususnya struktur
        beton.
            Sebagai informasi tambahan, dalam algoritma ini digunakan
        library Pandas. Sehingga pengguna harus memiliki mesin (komputer)
        yang sudah ter-install Pandas.
        
        Rizal Purnawan.
        
        """
        )
    
    def __init__(self, h, b, D, fc, fy, Mu, phi= 0.85, clear_cover= 40):
        # Keterangan:
        # h : tinggi balok beton (dalam mm)
        # b : lebar balok beton (dalam mm)
        # D : diameter tulangan (dalam mm)
        # fc : tegangan tekan ultimate beton (dalam MPa)
        # fy : tegangan leleh baja (dalam MPa)
        # Mu : momen terfaktor (dalam Nmm)
        # cc : bernilai 40 mm secara default,
        #      dan dapat diganti (dalam mm)
        self.h = h
        self.b = b
        self.D = D
        self.fc = fc
        self.fy = fy
        self.phi = phi
        self.Mu = Mu
        self.cc = clear_cover
        
        # Fungsi-fungsi bantuan:
        from math import sqrt, pi, ceil, floor
        self.__sqrt = sqrt
        self.__pi = pi
        self.__ceil = ceil
        self.__floor = floor
        
        # Atribut-atribut lain:
        self.Ac = h *b                  # luas penampang total balok.
        self.At = 0.25 *pi *D**2        # Luas penampang 1 batang tulangan.
        self.d = 0.8 *h                 # Jarak dari titik berat tulangan
                                        # tarik hingga serat ekstrim tertekan,
                                        # diasumsikan sebesar `0.8h` sebagai
                                        # permulaan.
        self.d_ = clear_cover + 0.5*D   # Jarak dari titik berat tulangan tekan
                                        # hingga serat ektrim tertekan.
        self.beta_1 = 0.85 if 17 <= fc <= 28 \
                        else ( 0.85 - 0.05*((fc - 28) / 7)
                              if 28 < fc < 55
                              else (0.65 if fc >= 55
                                    else None) 
                             )          # Nilai `beta_1` yang ditentukan
                                        # berdasarkan SNI 2847:2019
        self.As_min = ( max( [ sqrt(fc)/(4 *fy), 1.4/fy ] )
                        *b *self.d )    # Ini adalah nilai `As` minimum, yang
                                        # ditentukan berdasarkan SNI 2847:2019.
        self.As_max = ( 0.3825 *( fc/fy ) *self.beta_1
                        *b *self.d )    # Ini adalah nilai `As` maksimum.
        
    def __kand_r_dan_As(self, r_step= 0.05, As_step= 10):
        """
        Nama Metode    : __kand_r_dan_As
        Status         : Metode tersembunyi
        
        ___Tentang metode
        
            Metode ini berfungsi untuk menentukan nilai-nilai `r` dan
        `As` yang mungkin untuk digunakan. `r` adalah suatu rasio yang
        terdapat pada persamaan (8) dan (9) pada teori di atas.
        Nilai `r` ditentukan berdasarkan persamaan (11). Kemudian nilai
        `As` ditentukan berdasarkan persamaan (12). Output dari metode
        ini adalah instansiasi atribut dengan objek sebagaimana
        dijelaskan pada persamaan (13).
            Argumen `r_step` pada metode ini diberikan sejatinya bernilai
        0.05. Para pengguna dapat merubah nilai tersebut sesuai kebutuhan
        komputasi.
            Kemudian argumen `As_step` berfungsi serupa dengan `r_step`.
        Pada persamaan (12) memang tidak dijelaskan mengenai hal ini.
        Sehingga boleh saja himpunan `kand(As)` pada persamaan (12)
        berupa suatu kontinuum. Namun untuk keperluan komputasi, kita
        perlu melakukan diskritisasi, yaitu dengan menggunakan `As_step`.
            Sebagai tambahan, metode ini merupakan suatu metode
        tersembunyi, sehingga hanya dapat dipanggil di luar class dengan
        syntax:
        
            class._OptimasiLenturBalok__kand_r_As(parameter)
            
        ___Rizal Purnawan
        """
        kand_r = []
        r_0 = 0.2    # Hal ini karena nilai terkecil yg mungkin
                     # untuk `r` adalah 0.2.
        while True:
            # Kita lakukan iterasi berdasarkan persamaan (11) untuk
            # mendapatkan `kand_r`.
            r = r_0 + r_step
            if r > 0.6:
                break
                # Kita hentikan iterasi ketika nilai melewati 0.6,
                # sebagaimana dijelaskan pada persamaan (11)
            kand_r.append(r)
            r_0 = r
        kand_r_As = []    # List ini akan digunakan untuk
                          # menampung nilai output dari
                          # metode berdasarkan persamaan (13)
        Mu_per_phi = self.Mu / self.phi
        for r in kand_r:
            As_0 = self.As_min
            while True:
                As = As_0 + As_step
                if As > self.As_max:
                    break
                    # Kita hentikan iterasi ketika nilai melewati
                    # As_max.
                # Sekarang kita lakukan tes, apakah dengan nilai `As`
                # ini momen terpenuhi, sebagaimana dijelaskan pada
                # persamaan (12):
                M_tes = ( (self.d - r*self.d_) *self.fy *As
                          - (1 - r)**2 *( self.fy**2 *As**2
                                          )/( 1.7 *self.fc *self.b )
                          )    # Persamaan `M_tes` dihitung berdasarkan
                               # persamaan (9) pada teori di atas.
                if M_tes >= Mu_per_phi:
                    kand_r_As.append([r, As])
                As_0 = As
        self.kand_r_As = kand_r_As    # Kita instansiasikan `kand_r_As`
                                      # sebagai suatu atribut pada class.
            
    def __optimasi(self, n_opt= 20):
        """
        Nama Metode    : __optimasi
        Status         : Metode tersembunyi
        
        ___Tentang metode
        
            Metode ini berfungsi untuk melakukan optimasi terhadap
        nilai-nilai pada `self.kand_r_As` (yang telah dieksekusi pada
        metode __kand_r_As). Optimasi dilakukan dengan memilih `r` dan
        `As` terbaik sedemikian hingga nilai `As + As'` minimum. Untuk
        itu, kita perlu menghitung `a`, `fs` dan `As'`.
            Metode ini tidak memiliki argumen. Selain itu, metode ini
        juga merupakan suatu metode tersembunyi, sehingg hanya dapat
        diakses di luar class menggunakan syntax:
        
            class._OptimasiLenturBalok__optimasi()
            
        ___Rizal Purnawan
        """
        try:
            kand_r_As_a_Atot = []    # List ini disiapkan untuk nantinya
                                     # menampung nilai-nilai `r`, `As`, `a`
                                     # dan `As + As'`.
            for r, As in self.kand_r_As:
                a = ( (1 - r) *self.fy *As         # Nilai `a` berdasarkan
                     / (0.85 *self.fc *self.b) )        # persamaan (8).
                epsilon_s = ( 0.003 *self.d_       # Nilai `epsilon_s` ditetapkan
                              *self.beta_1 / a )   # sesuai persamaan (14), dengan
                                                   #`epsilon_cu = 0.3% = 0.003`.
                fs = ( epsilon_s *2*10**5          # Nilai `fs` sesuai
                       if epsilon_s < 0.002        # persamaan (15).
                       else self.fy )
                As_ = r * self.fy *As / fs         # Nilai `As'` ditetapkan sesuai
                                                   # persamaan (16).
                Atot = As + As_                    # Nilai `Atot` dihitung sebagai
                                                   # `As + As_`.
                kand_r_As_a_Atot.append([r, As, a, Atot])
            # Sekarang kita cari nilai `r`, `As` dan `a` pada `Atot` terkecil:
            try:
                opt_kand = sorted(kand_r_As_a_Atot, key= lambda i: i[-1])[:n_opt]
            except:
                opt_kand = sorted(kand_r_As_a_Atot, key= lambda i: i[-1])
            self.opt_kand = opt_kand
            self.best_r, self.best_As, self.best_a, self.best_As_ = opt_kand[0]
        except:
            print("ERROR: Eksekusi metode __kand_r_As terlebih dahulu!")
            
    def __delta(self, Nreb, Ntot, Nlay):
        """
        ___Tentang metode
        
            Metode ini merupakan metode bantuan, yang dibuat berdasarkan
        persamaan (21).
        """
        try:
            N_end = Ntot % Nreb
            if N_end == 0:
                N_k = [Nreb] *Nlay
            else:
                N_k = ( [Nreb] *(Nlay - 1) ) + [N_end]
            return (
                sum( [N_k[k - 1] *(self.cc
                                   + ( (k - 0.5) *self.D )
                                   + ( 25 *k )
                                   - 25)
                      for k in range(1, Nlay + 1)]
                   ) / Ntot,
                N_k )
        except:
            print("ERROR: Argumen TIDAK valid!")
            
    def __evaluasi(self, cari= False):
        """
        Nama Metode    : __evaluasi
        Statu          : Metode tersembunyi
        
        ___Tentang metode
        
            Metode ini digunakan untuk melakukan evaluas terhadap
        atribut-atribut hasil optimasi (setelah metode __optimasi
        dieksekusi). Tujuan utama dari metode ini adalah memastikan
        bahwa dengan `As` dan `As'` hasil optimasi, momen nominal
        aktual `Mn` memenuhi, yaitu
        
            phi Mn >= Mu .
            
            Argumen `cari` pada metode ini dapat ditentukan sebagai
        `True` atau `False`. Jika bernilai `True`, maka akan dilakukan
        iterasi pada `self.opt_kand` untuk mencari ulang nilai 'As'
        dan `As'` terbaik sedemikian hingga kontrol momen nominal
        aktual terpenuhi. Sejatinya, kami tetapkan `cari` bernilai
        `False`.
            Selain itu, metode ini juga akan memunculkan konfigurasi
        tulangan tarik maupun tekan.
        
        ___Rizal Purnawan
        """
        try:
            N_reb = self.__floor(              # Nilai `N_reb` ditentukan
                ( self.b + 25 - 2 *self.cc )   # berdasarkan
                / ( 25 + self.D )              # persamaan (17).
            )
            if cari == True:
                opt_kand = self.opt_kand
            elif cari == False:
                opt_kand = [self.opt_kand[0]]
            for r, As, _, Atot in opt_kand:
                As_ = Atot - As    # Nilai `As'`.
                Ntot, Nlay, Nend = list(), list(), list()
                for A in [As, As_]:
                    _Ntot = self.__ceil(             # Menentukan `Ntot`
                        ( 4 *A )                     # berdasarkan
                        / ( self.__pi *self.D**2 )   # persamaan (18)
                    )
                    _Ntot = max([_Ntot, 2])
                    Ntot.append(_Ntot)
                    _Nlay = self.__ceil(             # Menentukan 'Nlay'
                        _Ntot / N_reb                # berdasarkan
                    )                                # persamaan (19)
                    Nlay.append(_Nlay)
                    _Nend = _Ntot % N_reb    # `Nend` ditentukan berdasarkan
                                             # persamaan (20)
                    Nend.append(_Nend)
                # Kita hitung ulang `As` dan `As_`:
                As = Ntot[0] *0.25 *self.__pi *self.D**2
                As_ = Ntot[1] *0.25 *self.__pi *self.D**2
                # Sekarang kita tentukan nilai `d` berdasarkan persamaan (22)
                delta_tar = self.__delta(N_reb, Ntot[0], Nlay[0])
                d = self.h - delta_tar[0]
                
                # Dan nilai `d_` ditentukan berdasarkan persamaan (23)
                delta_tek = self.__delta(N_reb, Ntot[1], Nlay[1])
                d_ = delta_tek[0]
                
                fs = r *self.fy *As / As_    # Ditentukan `fs` berdasarkan
                                             # `As` dan `As_` yg telah direvisi
                a = (self.fy *As - fs *As_) / (0.85 *self.fc *self.b)
                Mn = ( (self.fy *As *(d - 0.5*a))    # `Mn` dihitung berdasarkan
                       + (fs *As_ *(0.5*a - d_)) )   # persamaan (1)
                if self.phi *Mn >= self.Mu:
                    # Jika persamaan (1) telah terpenuhi, maka kita simpan
                    # informas-informasi yang kita dapatkan dan kita
                    # akhiri iterasi:
                    self.tulangan_tarik = f"{Ntot[0]}D{int(self.D)}"
                    self.tulangan_tekan = f"{Ntot[1]}D{int(self.D)}"
                    self.konf_tarik = delta_tar[1]
                    self.konf_tekan = delta_tek[1]
                    self.Mn = Mn
                    self.As = As
                    self.As_ = As_
                    break
        except:
            print("ERROR: Eksekusi metode __optimasi terlebih dahulu!")
            
    def proses(self, r_step= 0.05, As_step= 10, n_opt= 20, cari= False):
        """
        Nama Metode    : proses
        Status         : Dapat diakses
        
        ___Tentang metode
        
            Metode ini merupakan unifikasi dari seluruh metode yang
        telah dibuat. Metode ini dapat diakses dari luar class. Untuk
        mengetahui metode-metode yang dieksekusi oleh metode ini
        serta mengenai informasi metode-metode tersebut, jalankanlah
        program berikut:
        
            class._OptimasiLenturBalok__kand_r_dan_As.__doc__
            class._OptimasiLenturBalok__optimasi.__doc__
            class._OptimasiLenturBalok__evaluasi.__doc__
            
        ___Rizal Purnawan
        """
        try:
            import pandas as pd
            self.__kand_r_dan_As(r_step= r_step, As_step= As_step)
            self.__optimasi(n_opt= n_opt)
            self.__evaluasi(cari= cari)
            ket = "AMAN" if self.phi * self.Mn >= self.Mu else "TIDAK AMAN!"
            output_dic = {"Tulangan Tarik": [self.tulangan_tarik, " "],
                          "Tulangan Tekan": [self.tulangan_tekan, " "],
                          "Konfigurasi Tulangan Tarik":
                              [self.konf_tarik, f"{len(self.konf_tarik)} Baris"],
                          "Konfigurasi Tulangan Tekan":
                              [self.konf_tekan, f"{len(self.konf_tekan)} Baris"],
                          "Luas Tulangan Tarik": [self.As, "(mm-psg)"],
                          "Luas Tulangan Tekan": [self.As_, "(mm-psg)"],
                          "Momen Nominal Terfaktor":
                              [self.phi *self.Mn / 10**6, "(kNm)"],
                          "Momen Ultimate": [self.Mu / 10**6, "(kNm)"],
                          "Keterangan": [ket, " "]}
            df = pd.DataFrame(output_dic)
            output_df = df.T
            output_df.columns = ["Nilai", "Satuan"]
            display(output_df)
        except:
            print("ERROR: Input TIDAK VALID!, atau Penampang TIDAK memenuhi!")
