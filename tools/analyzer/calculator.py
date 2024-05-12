import numpy as np
from analyzer.anphonio import ParseResult

class Calculator:
    def __init__(self, file_result_3ph,
                 file_result_4ph=None,
                 file_isotope=None,
                 average_gamma=True):
        self.file_result_3ph = file_result_3ph
        self.file_result_4ph = file_result_4ph
        self.file_isotope = file_isotope
        self.omega = None  # Frequency array
        self.omega4 = None
        self.gamma3 = None    # linewidth due to 3-phonon scattering
        self.gamma4 = None    # linewidth due to 4-phonon scattering
        self.gamma4_interpolated = None
        self.gamma_iso = None  # linediwth due to isotope scattering
        self.vel = None    # Velocity array
        self.vel4 = None
        self.qpoint_weight = None  # Weight array
        self.qpoint_weight4 = None
        self.qpoints = None
        self.qpoints4 = None
        self.qgrid = None
        self.qgrid4 = None
        self.volume = None
        self.volume4 = None
        self.temperatures = None
        self.temperatures4 = None
        self.average_gamma = average_gamma

        self._BOHR = 0.52917721092
        self._k_Boltzmann = 1.3806488e-23

        if self.file_result_3ph is not None:
            self.set_variables_3ph()

        if self.file_result_4ph is not None:
            self.set_variables_4ph()

        if self.file_isotope is not None:
            self.set_variables_iso()

    def set_variables_3ph(self):
        result = ParseResult(self.file_result_3ph)
        self.omega = result.omega
        if self.average_gamma:
            self.gamma3 = self.average_gamma_at_degenerate_point(self.omega, result.gamma)
        else:
            self.gamma3 = result.gamma
        self.vel = result.vel
        self.volume = result.volume
        self.qpoint_weight = result.multiplicity
        self.qpoints = result.q_coord
        self.qgrid = result.kgrid
        self.temperatures = result.temperatures

    def set_variables_4ph(self):
        result = ParseResult(self.file_result_4ph)
        self.omega4 = result.omega
        if self.average_gamma:
            self.gamma4 = self.average_gamma_at_degenerate_point(self.omega4, result.gamma)
        else:
            self.gamma4 = result.gamma
        self.vel4 = result.vel
        self.volume4 = result.volume
        self.qpoint_weight4 = result.multiplicity
        self.qpoints4 = result.q_coord
        self.qgrid4 = result.kgrid
        self.temperatures4 = result.temperatures

    def set_variables_iso(self):
        result = np.loadtxt(self.file_isotope)
        nk = round(np.max(result[:,0]))
        ns = round(np.max(result[:,1]))
        omega = result[:,2].reshape((nk, ns))
        gamma = result[:,3].reshape((nk, ns, 1))
        if self.average_gamma:
            self.gamma_iso = self.average_gamma_at_degenerate_point(omega, gamma)[:,:,0]
        else:
            self.gamma_iso = gamma[:,:,0]

    def check_data_load(self, four_phonon, isotope):
        if self.gamma3 is None:
            self.set_variables_3ph()

        if four_phonon and (self.gamma4 is None):
            if self.file_result_4ph is None:
                raise RuntimeError("file_result_4ph must be given when initializing the class")
            else:
                self.set_variables_4ph()

        if isotope and (self.gamma_iso is None):
            if self.file_isotope is None:
                raise RuntimeError("file_isotope must be given when initializing the class")
            else:
                self.set_variables_iso()

    def average_gamma_at_degenerate_point(self, frequencies, gamma, tol_omega=1e-3):
        """
        Averages the gamma values at degenerate points across the phonon band structure.
        :param tol_omega: Tolerance for considering frequencies as degenerate.
        """
        nk, ns = frequencies.shape
        nt = gamma.shape[-1]

        gamma_avg = np.zeros((nk, ns, nt), dtype=float)

        # Loop over all k-points
        for i in range(nk):
            degeneracy_at_k = []
            omega_prev = frequencies[i, 0]
            ideg = 1

            # Identify degenerate modes
            for j in range(1, ns):
                omega_now = frequencies[i, j]
                if abs(omega_now - omega_prev) < tol_omega:
                    ideg += 1
                else:
                    degeneracy_at_k.append(ideg)
                    ideg = 1
                    omega_prev = omega_now
            degeneracy_at_k.append(ideg)  # Append the last set of degeneracies

            is_index = 0

            # Average the tau values for each set of degenerate modes
            for deg in degeneracy_at_k:
                damp_sum = np.zeros(nt)
                for k in range(is_index, is_index + deg):
                    damp_sum += gamma[i, k, :]

                for k in range(is_index, is_index + deg):
                    gamma_avg[i, k, :] = damp_sum / deg

                is_index += deg

        return gamma_avg

    def compute_linewidth(self, temperature, four_phonon=False, isotope=False):

        self.check_data_load(four_phonon, isotope)

        tempdiff = np.abs(self.temperatures - temperature)
        index = np.argsort(tempdiff)

        if tempdiff[index[0]] > 0:
            raise RuntimeWarning("The data at exactly at {} K not found. "
                                 "Return the values at {} K instead".format(temperature,
                                                                            self.temperatures[index[0]]))

        return self.gamma3[:,:,index[0]]

    def calc_kappa(self):
        # Placeholder for translating the calc_kappa function
        pass

    def heat_capacity(self, omega, temp):
        return omega / temp

    # Additional methods will be defined here.
