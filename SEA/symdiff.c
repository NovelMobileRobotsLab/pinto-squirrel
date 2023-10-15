
#include <math.h>
#include <stdio.h>


int main(int argc, char const *argv[]) {
    // x0 = pow(T2, 2);
    // x1 = pow(tau_max, 2);
    // x2 = 4 * x1;
    // x3 = pow(omega_max, 2);
    // x4 = m_l * x3;
    // x5 = pow(T1, 4) * x0;
    // x6 = x2 + x4 * x5;
    // x7 = sqrt(m_l * x6);
    // x8 = pow(T1, 2);
    // x9 = omega_max * x8;
    // x10 = T2 * x9;
    // x11 = x10 * x7;
    // x12 = m_l * x10;
    // x13 = m_l * x2;
    // x14 = pow(m_l, 2);
    // x15 = x3 * x5;
    // x16 = sqrt(x13 + x14 * x15);
    // x17 = x12 + x16;
    // x18 = 2 * tau_max;
    // x19 = x18 / x17;
    // x20 = 1 / (T2 * m_l);
    // x21 = (1.0 / 2.0) * t / tau_max;
    // x22 = x20 * x21;
    // x23 = x4 * x8 * exp(x17 * x22);
    // x24 = x20 * x7;
    // x25 = T2 * m_l;
    // x26 = T2 / x7;
    // x27 = x12 - x16;
    // x28 = x18 / x27;
    // x29 = x8 * exp(x22 * x27);
    // x30 = pow(tau_max, 4);
    // x31 = pow(T1, 6) * pow(T2, 3) * pow(omega_max, 3);
    // x32 = 2 * x30;
    // x33 = pow(T1, 8) * pow(T2, 4) * pow(omega_max, 4) * x14;
    // x34 = x13 * x15;
    // x_acode = omega_max * t + x0 * x14 * x28 * x29 * x3 * (-m_l * x31 * x7 - 2 * x1 * x11 + x32 + x33 + x34) / (pow(T1, 10) * pow(T2, 5) * pow(m_l, 3) * pow(omega_max, 5) + 6 * x1 * x14 * x31 + 8 * x12 * x30 - x32 * x7 - x33 * x7 - x34 * x7) - x0 * x19 * x23 * (x12 - x7) / (-x11 + x6) + x23 * x26 * ((x17 != 0) ? (-x19 * x25 * exp(-x21 * (x24 + x9))) : (t)) - x26 * x29 * x4 * ((x27 != 0) ? (-x25 * x28 * exp(-x21 * (-x24 + x9))) : (t));

    float x = 0.1;
    for(int i = 0; i < 100; i++){
        x = pow(0.1, 2);
    }

    printf("%f", x);
    return 0;
}
