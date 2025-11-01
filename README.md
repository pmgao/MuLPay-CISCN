# MuLPay-CISCN
MuLPay：基于格的高效后量子安全区块链钱包系统 (第十六届全国大学生信息安全竞赛)

## 编译与执行方法

### 后端部分

```shell
$ cd ./algorithm

$ nvcc -O3 -arch=sm_XX SM4_CTR.cu -o encrypt.exe # XX is determined by the model of GPU

$ nvcc -O3 -arch=sm_XX SM4_CTR.cu -o decrypt.exe # XX is determined by the model of GPU

$ mv encrypt.exe decrypt.exe ../app/

$ cd ../app/

$ sage app.sage
```

### 前端部分

``` 
http://127.0.0.1:5000/ 
```

## 技术细节	

MuLPay具备社交恢复、共同支付和单人管理等功能，且提供了完整的前端用户交互界面。对应报告PDF可见：[Report](https://github.com/pmgao/MuLPay-CISCN/blob/main/Report.pdf).

参与者：高培民（山东大学，香港理工大学），吴兵（山东大学，香港理工大学）、欧阳仁鼎（山东大学，浙江大学）、张雨欣（山东大学，上海交通大学）。

---------------------------------------------------------------------------------------------------------------------

MuLPay: Efficient Lattice-Based Post-Quantum Secure Blockchain Wallet System (The 16th National College Student Information Security Competition)
## Compile and execution

### Back End

```shell
$ cd ./algorithm

$ nvcc -O3 -arch=sm_XX SM4_CTR.cu -o encrypt.exe # XX is determined by the model of GPU

$ nvcc -O3 -arch=sm_XX SM4_CTR.cu -o decrypt.exe # XX is determined by the model of GPU

$ mv encrypt.exe decrypt.exe ../app/

$ cd ../app/

$ sage app.sage
```

### Front End

``` 
http://127.0.0.1:5000/ 
```

## Technical Details

MuLPay has functions such as social recovery, joint payment, and single-person management, and provides a complete front-end user interface. The corresponding report PDF can be found at [Report](https://github.com/pmgao/MuLPay-CISCN/blob/main/Report.pdf).

Members: Peimin Gao (Shandong University, The HK PolyU), Bing Wu (Shandong University, The HK PolyU), Rending Ouyang (Shandong University, ZJU), and Yuxin Zhang (Shandong University, SJTU).

