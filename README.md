# Wifi-System-Performance-Analysis

IEEE 802.11, commonly known as WiFi, is a famous wireless communication protocol. In this final project, we aim to simulate WiFi and compare the system performance under different scenarios. First, We will evaluate the system performance of standalone WiFi. Then, we will evaluate the system performance of WiFi when 4G LTE coexists. Finally, we will evaluate the system performance of a hybrid system with uplink following WiFi and downlink following 4G. We find out that different architectures are suitable for different scenarios.

**Please refer to report.pdf for detailed information.**
# How to execute
There are 3 subdirectories:
- hybrid/
- coexistence/
- wifi/
 
In /hybrid/, there are codes for simulation in section 4.3 (Hybrid system).
execute
  - hybrid.m: scenario
  - hybrid_different_lambda.m: comparison of different value of lambda value
  - hybrid_different_ms.m: comparison of different number of MS per cell

In /src_code/coexistence/, there are codes for simulation in section 4.2(Coexistence of WiFi and 4G LTE).
execute
- WiFi_LTE.m: scenario
- WiFi_LTE_with_diff_lambda.m: comparison of different value of lambda value
- WiFi_LTE_with_diff_ms.m: comparison of different number of MS per cell

In /src_code/wifi/, there are codes for simulation in section 4.1 (Standalone WiFi).
execute
- WiFi.m: scenario
- WiFi_different_lambda.m: comparison of different value of lambda value
- WiFi_different_ms.m: comparison of different number of MS per cell
- WiFi_different_CWmin.m: comparison of different value of CWmin
- WiFi_different_Ith.m: comparison of different interference threshold

