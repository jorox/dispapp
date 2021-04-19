# Post-processing

## The log files

The log files produced by fix-phonon are not very hard to understand. A quick look at src/fix-phonon.cpp shows that:


    // write log file, here however, it is the dynamical matrix that is written
    fprintf(flog,"############################################################\n");
    fprintf(flog,"# Current time step                      : " BIGINT_FORMAT "\n", update->ntimestep);
    fprintf(flog,"# Total number of measurements           : %d\n", neval);
    fprintf(flog,"# Average temperature of the measurement : %lg\n", TempAve);
    fprintf(flog,"# Boltzmann constant under current units : %lg\n", boltz);
    fprintf(flog,"# basis vector A1 = [%lg %lg %lg]\n", basevec[0], basevec[1], basevec[2]);
    fprintf(flog,"# basis vector A2 = [%lg %lg %lg]\n", basevec[3], basevec[4], basevec[5]);
    fprintf(flog,"# basis vector A3 = [%lg %lg %lg]\n", basevec[6], basevec[7], basevec[8]);
    fprintf(flog,"############################################################\n");
    fprintf(flog,"# qx\t qy \t qz \t\t Phi(q)\n");

    EnforceASR();

    // to get D = 1/M x Phi
    for (idq = 0; idq < ntotal; ++idq) {
      ndim =0;
      for (idim = 0; idim < fft_dim; ++idim)
      for (jdim = 0; jdim < fft_dim; ++jdim) Phi_all[idq][ndim++] *= M_inv_sqrt[idim/sysdim]*M_inv_sqrt[jdim/sysdim];
    }

    idq =0;
    for (int ix = 0; ix < nx; ++ix) {
      double qx = double(ix)/double(nx);
      for (int iy = 0; iy < ny; ++iy) {
        double qy = double(iy)/double(ny);
        for (int iz = 0; iz < nz; ++iz) {
          double qz = double(iz)/double(nz);
          fprintf(flog,"%lg %lg %lg", qx, qy, qz);
          for (idim = 0; idim < fft_dim2; ++idim)
            fprintf(flog, " %lg %lg", std::real(Phi_all[idq][idim]),
                                      std::imag(Phi_all[idq][idim]));
          fprintf(flog, "\n");
          ++idq;
        }
      }
    }
    fflush(flog);    

And sample generated log file


    # Current time step                      : 6500000
    # Total number of measurements           : 600000
    # Average temperature of the measurement : 299.371
    # Boltzmann constant under current units : 8.61734e-05
    # basis vector A1 = [2.56873 0 0]
    # basis vector A2 = [1.28436 2.22458 0]
    # basis vector A3 = [1.28436 0.741528 2.09736]
 
    # qx	 qy 	 qz 		 Phi(q)
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0 0 -0
    0 0 0.125 0.00585736 6.8795e-25 -0.000156973 -9.29854e-05 -7.36102e-05 1.46466e-06 -0.000156973 9.29854e-05 0.0055872 0 -3.49646e-05 -1.0267e-06 -7.36102e-05 -1.46466e-06 -3.49646e-05 1.0267e-06 0.0341953 1.48916e-25
    0 0 0.25 0.0190152 0 7.96227e-05 0.000201581 -0.000219957 7.29924e-06 7.96227e-05 -0.000201581 0.0187525 0 -0.000154354 1.0039e-05 -0.000219957 -7.29924e-06 -0.000154354 -1.0039e-05 0.112138 -1.28116e-23
    0 0 0.375 0.0313595 -2.7518e-24 0.000281555 -0.000344008 -0.000342451 -8.86663e-06 0.000281555 0.000344008 0.0320947 -2.74977e-24 -0.000221315 -1.10892e-05 -0.000342451 8.86663e-06 -0.000221315 1.10892e-05 0.186243 1.76698e-23

The code shows that what is being written to the log file is the Φ function itself (unlike for the binary file actually!). Great, less work for us!

If we look at an example log file we can see that the first three values are the coordinates of q and the rest are the real and complex parts for the terms of the Φ matrix.

So for a primitive unit cell, for example, we can expect 3x3x2=19 values.
No need to worry about whether it's column or row sorted since the eigenvalues of the transpose and the matrix are the same.

## phana.ipynb

The generated data is post-processed using a custom Python module which reads the log files produced by fix-phonon.
The module returns the last dynamical matrix which seems to be the same thing that phana does.
Then the dispersion curve is calculated along a selected path of special points.
The data is interpolated using cubic splines and the eigenvalues are sorted.

## Plotting

For plotting an interactive figure we use Bokeh with custom JS callback functions