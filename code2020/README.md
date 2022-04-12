# Xray_astronomy
Course assignment
</p>
If any questions, please contact me.(Wechat or my mailbox(baotong@smail.nju.edu.cn))

# Analysis 
Read "xmm_abc_guide.pdf" first(chapter 5&6).

Run "run_XMMproducts_spectra.py" to process the data and get spectra.

Run "xmm_fits_to_txt.py" and "pfold_xmm.py" in order to get phase-folded light curve.

感谢孙逸辰同学的提醒，由于信噪比非常高，背景level在光变曲线中太淡了看不出来。建议将pfold_xmm.py第114行中的颜色改为blue,alpha改为0.9。

Note: 请下载最新的code，之前的pfold_xmm.py中有一些偏差，会影响light curve的精准度。
# Fit the spectra with model
Query "XspecManual.pdf" for help.

Fit the spectra with mkcf with absorption.

Read xspec_abc.txt for detailed.

# 作业提交
提交的光谱信息为图像（**.ps），参数拟合结果（截个图就行）以及highT的误差范围（示例如（4.23，10.78))

另外请附上你的光谱文件以及响应文件。(*.pha, *.arf, *.rmf)

请将上述文件做成压缩包后，发送到我的邮箱：baotong@smail.nju.edu.cn





