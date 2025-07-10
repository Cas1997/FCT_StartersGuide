std::function<void(const double*, double*)> field() {
  return [](const double* x, double* b) {
    if( x[2] > 440 && x[2] < 500 && (x[0] * x[0] + x[1] * x[1]) < 50*50){
        b[0] = 0.; // In kGauss. 1 kGauss = 0.1 T
        b[1] = 2.5;
        b[2] = 0.;
    } else {
        b[0] = 0.; // In kGauss
        b[1] = 0.;
        b[2] = 20.;
    }
  };
}