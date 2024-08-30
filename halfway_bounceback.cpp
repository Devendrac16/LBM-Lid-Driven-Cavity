///////Devendra Chaudhari 
//////LId Driven Cavity 
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
//relaxation time calculation
double Relaxation_time(double Re, double Ma, int L)
{
    double U, kin_visc, relax_time;
    double cs = 1 / sqrt(3.0);
    U = Ma * cs;
    kin_visc = U * L / Re;
    relax_time = 3 * kin_visc + 0.5;
    return relax_time;
}

int main()
{
    double Re = 400, Ma = 0.5, tau, ux, uy, rho, u, term1, term2;
    int Nx, Ny, ix, iy;
    const int q = 9;
    const double rho0 = 1.0;
    double cs = 1 / sqrt(3.0);
    int iteration = 0;
    double error = 0;
    // double uw = 1.0; 
    double uw =  Ma*cs;
    // double U0= Ma*1 / sqrt(3.0);
    std::cout<<"enter number of nodes in x-direction: ";
    std::cin>>Nx;
    std::cout<<"enter number of nodes in y-direction: ";
    std::cin>>Ny;

    std::vector<std::vector<double>> ux_old(Nx + 3, std::vector<double>(Ny + 3, 1.0));
    std::vector<std::vector<double>> uy_old(Nx + 3, std::vector<double>(Ny + 3, 1.0));
    std::vector<std::vector<double>> ux_final(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    std::vector<std::vector<double>> uy_final(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    std::vector<std::vector<double>> u_final(Nx + 3, std::vector<double>(Ny + 3, 0.0));
    std::vector<std::vector<double>> rho_final(Nx + 3, std::vector<double>(Ny + 3, rho0));
    std::vector<std::vector<double>> p(Nx+3,std::vector<double>(Ny+3,0.0));

    std::setprecision(7);

    tau = Relaxation_time(Re, Ma, (Nx-1));

    int cx[q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    double w[q] = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};
    int opposite[q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    std::vector<std::vector<std::vector<double>>> f(q, std::vector<std::vector<double>>(Nx + 3, std::vector<double>(Ny + 3, 0.0)));
    std::vector<std::vector<std::vector<double>>> fcoll(q, std::vector<std::vector<double>>(Nx + 3, std::vector<double>(Ny + 3, 0.0)));
    std::vector<double> feq(q, 0.0);

    // Initialization polpulation

    for (int i = 1; i <= Nx+1; i++)
    {
        for (int j = 1; j <= Ny+1; j++)
        {
            for (int a = 0; a < q; ++a)
            {
                f[a][i][j] = w[a] * rho0;
            }
        }
    }

    std::ofstream file;
    file.open("lbm_output.dat");

    std::ofstream data;
    data.open("error.txt");
    data<<"x-grid ponts"<<Nx+1<<"\t"<<"Y-grid points"<<Ny+1<<"\n";
    data<<"i"<<"j"<<"\t\t\t"<<"rho\t\t\t"<<"ux\t\t\t"<<"uy\t\t\t"<<"u\t\t\t"<<"psi\t\t\t"<<"omega\t\t\t"<<"error\n";

    std::ofstream ux_vel;
    ux_vel.open("ux_velocity.dat");
    

    std::ofstream uy_vel;
    uy_vel.open("uy_velocity.dat");
   

    do
    {
        //////////////////collision///////////////////////////////////////
        for (int i = 1; i <= Nx+1; i++)
        {
            for (int j = 1; j <= Ny+1; j++)
            {
                ux = ux_final[i][j];
                uy = uy_final[i][j];
                rho = rho_final[i][j];
                for (int a = 0; a < q; a++)
                {
                    term1 = ux * cx[a] + uy * cy[a];
                    term2 = term1 * term1;
                    feq[a] = w[a] * rho * (1.0 + 3.0 * term1 + 4.5 * term2 - 1.5 * (ux * ux + uy * uy));
                    fcoll[a][i][j] = f[a][i][j] - (f[a][i][j] - feq[a]) / tau; // collision
                }
            }
        }
        
        ////////////////////////streaming/////////////////////////
        for (int i = 1; i <= Nx+1; i++)
        {
            for (int j = 1; j <= Ny+1; j++)
            {
                for (int a = 0; a < q; a++)
                {
                    int nextX = i + cx[a];
                    int nextY = j + cy[a];
                    f[a][nextX][nextY] = fcoll[a][i][j]; 
                }
            }
        }
        //Boundary conditions: Top lid moves with constant velocity uw=1.0
        for (int i = 1; i <= Nx+1; i++)
        {
            f[4][i][Ny+1] = fcoll[opposite[4]][i][Ny+1];
            f[7][i][Ny+1] = fcoll[opposite[7]][i][Ny+1] - 6.*w[opposite[7]]*rho0*uw ;
            f[8][i][Ny+1] = fcoll[opposite[8]][i][Ny+1] + 6.*w[opposite[8]]*rho0*uw ;

            f[2][i][1] = fcoll[opposite[2]][i][1];
            f[5][i][1] = fcoll[opposite[5]][i][1];
            f[6][i][1] = fcoll[opposite[6]][i][1];
        }

        for (int j = 1; j <= Ny+1; j++)
        {
            f[1][1][j] = fcoll[opposite[1]][1][j];
            f[5][1][j] = fcoll[opposite[5]][1][j];
            f[8][1][j] = fcoll[opposite[8]][1][j];

            f[3][Nx+1][j] = fcoll[opposite[3]][Nx+1][j];
            f[6][Nx+1][j] = fcoll[opposite[6]][Nx+1][j];
            f[7][Nx+1][j] = fcoll[opposite[7]][Nx+1][j];
        }


        // Corner Treatment
        // Bottom-left corner (1,1)
        f[1][1][1] = fcoll[opposite[1]][1][1];
        f[2][1][1] = fcoll[opposite[2]][1][1];
        f[5][1][1] = fcoll[opposite[5]][1][1];

        // Bottom-right corner (Nx+1,1)
        f[3][Nx+1][1] = fcoll[opposite[3]][Nx+1][1];
        f[2][Nx+1][1] = fcoll[opposite[2]][Nx+1][1];
        f[6][Nx+1][1] = fcoll[opposite[6]][Nx+1][1];

        // Top-left corner (1,Ny+1)
        f[1][1][Ny+1] = fcoll[opposite[1]][1][Ny+1];
        f[4][1][Ny+1] = fcoll[opposite[4]][1][Ny+1];
        f[8][1][Ny+1] = fcoll[opposite[8]][1][Ny+1]+6.*w[opposite[8]]*rho0*uw;

        // Top-right corner (Nx+1,Ny+1)
        f[3][Nx+1][Ny+1] = fcoll[opposite[3]][Nx+1][Ny+1];
        f[4][Nx+1][Ny+1] = fcoll[opposite[4]][Nx+1][Ny+1];
        f[7][Nx+1][Ny+1] = fcoll[opposite[7]][Nx+1][Ny+1]- 6.*w[opposite[7]]*rho0*uw;


        /////////////////////macros////////////////
        double sum_unew = 0.;
        double sum_uold = 0.;
        for (int i = 1; i <= Nx+1; i++)
        {
            for (int j = 1; j <= Ny+1; j++)
            {
                ux = 0.0;
                uy = 0.0;
                rho = 0.0;
                for (int a = 0; a < q; a++)
                {
                    rho += f[a][i][j];
                    ux += f[a][i][j] * cx[a];
                    uy += f[a][i][j] * cy[a];
                }
                ux /= rho;
                uy /= rho;
                ux_final[i][j] = ux;
                uy_final[i][j] = uy;
                double t= (ux*ux)+(uy*uy);
                u_final[i][j]= sqrt(t);
                rho_final[i][j] = rho;
                p[i][j] = ((rho) * (u_final[i][j]));

                sum_unew += (ux - ux_old[i][j])*(ux - ux_old[i][j]) + (uy - uy_old[i][j])*(uy - uy_old[i][j]);
                sum_uold += ux_old[i][j]*ux_old[i][j] + uy_old[i][j]*uy_old[i][j];

                ux_old[i][j] = ux;
                uy_old[i][j] = uy;
            }
        }
        error = sqrt(sum_unew/sum_uold);

        if(iteration % 100 == 1)
        {
            std::cout << "Iteration " << iteration << ": Error = " << std::setprecision(7) << error << "\n";
        }

        iteration++;
    
    } while (error > 1e-6);
    std::vector<std::vector<double>>psi(Nx+3, std::vector<double>(Ny+3,0.0));

    for (int i = 1; i <= Nx+1; i++)
    {
        psi[i][0] = psi[i - 1][0] + 1.0 * uy_final[i][0];
    }
    for (int j = 1; j <= Ny+1; j++)
    {
        psi[0][j] = psi[0][j - 1] - 1.0 * ux_final[0][j];
    }
    for (int i = 1; i < Nx+1; i++)
    {
        for (int j = 1; j < Ny+1; j++)
        {
            psi[i][j] = psi[i - 1][j] + 1.0 * uy_final[i][j];
        }
    }

    std::vector<std::vector<double>>omega(Nx+3,std::vector<double>(Ny+3,0.0));
    // Calculate vorticity
        for (int i = 1; i < Nx; i++) {
            for (int j = 1; j < Ny; j++) {
                double dudy = (ux_final[i][j + 1] - ux_final[i][j - 1]) / (2.0 * 1.0);
                double dvdx = (uy_final[i + 1][j] - uy_final[i - 1][j]) / (2.0 * 1.0);
                omega[i][j] = dvdx - dudy;
            }
        }
        for(int i=1;i<=Nx+1;i++)
            {
                for(int j=0;j<=Ny+1;j++)
                {
                    data<<i<<j<<"\t\t"<<rho_final[i][j]<<"\t\t"<<ux_final[i][j]<<"\t\t"<<uy_final[i][j]<<"\t\t"<<u_final[i][j]<<"\t\t"<<psi[i][j]<<"\t\t"<<omega[i][j]<<"\n";
                }
            }

        file << "TITLE = \"LBM Lid-Driven Cavity\"\n";
        file << "VARIABLES = \"X\", \"Y\", \"Ux\", \"Uy\", \"Rho\",\"psi\",\"omega\", \"p\"\n";
        file << "ZONE I=" << Nx + 1 << ", J=" << Ny + 1 << ", F=POINT\n";

        for (int i = 1; i <= Nx+1; i++)
        {
            for (int j = 1; j <= Ny+1; j++)
            {
                file << i << " " << j << " " << ux_final[i][j] << " " << uy_final[i][j] << " " << rho_final[i][j]*(1./3.) << " " << psi[i][j] << " " << omega[i][j] << " " << p[i][j] << " " << "\n";
            }
        }
        for(int i=1; i<= Nx+1;i++)
        {
            double i1=(double)  Nx;
            double xaxis= (i-1)/i1;
            uy_vel<<xaxis<<"\t"<<(uy_final[i][(Nx-1)/2]/uw)<<"\n"; 
        }
        for(int j = 1; j <= Ny+1; j++)
        {
            double i2= (double) Ny;
            double yaxis= (j-1)/ i2;
            ux_vel<<"\t"<<(ux_final[(Ny-1)/2][j]/uw)<<"\t"<<yaxis<<"\n";
        }

    ux_vel.close();
    uy_vel.close();
    file.close();
    data.close();
    return 0;
}
