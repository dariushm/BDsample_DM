
double calc_psi(Fullmol *bases);
double calc_psi_2pi(Fullmol *bases, double leglen, int it);
double calc_psi_2pi(Fullmol *bases, double leglen1, double leglen0);
double calc_psi_1pi(Fullmol *bases, double leglen1, double leglen0);
void calc_three_angle(double &cthet1, double &cthet2, double &dih, Fullmol *bases, double dx, double dy, double dz, double R1, double leglen1, double leglen0, double *v, double *v1, double *n1, double *n2);
double numer_calc_Jacobian(Fullmol *bases, double psi0,  double leglen1, double leglen0, int molrot);

double set_norm_to_psi(double psi0, Fullmol *bases, double leglen1, double leglen0, int molfix, int molrot, double *csettoim);

double rotate_phi2_mol(double psi0, Fullmol *bases, double leglen1, double leglen0, int molfix, int molrot, double *csettoim);
