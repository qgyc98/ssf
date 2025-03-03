#ifndef SLIPSURF_H
#define SLIPSURF_H

long   detect_plastic_zones(long *plast_ip, double gamma_min);
long   detect_plast_ip(long nadjelem, long *adjelem, long *plast_ip, long zone_id, long *plast_adjelem);
double detect_max_gamma(long *plast_ip, long &max_gamma_ip, long &max_gamma_eid, double gamma_min);
long   approx_slip_surf(long id, long *plast_ip, double &a, double &b, double &minx, double &maxx);


#endif
