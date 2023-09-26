This package provides  test demo of our implementation for paper:
"VITON-DR: Details Retention Virtual Try-on via Non-rigid Registration "    

##Installation
`python 3.8.11`  
`pytorch 1.9,0`  
`numpy 1.21.5`
`opencv-python 4.6.0`  
`MATLAB 2021b`  


##test     
`python test.py`  

##Checkpoint 
We used the pre-trained models `latest_net_G.pth`, `latest_net_G1.pth`, `latest_net_G2.pth` from ACGPN, which can be download [here](https://drive.google.com/file/d/1UWT6esQIU_d4tUm8cjxDKMhB8joQbrFx/view?usp=sharing) .

## Dataset
**VITON Dataset** This dataset is presented in [VITON](https://github.com/xthan/VITON), containing 19,000 image pairs, each of which includes a front-view woman image and a top clothing image. After removing the invalid image pairs, it yields 16,253 pairs, further splitting into a training set of 14,221 pairs and a testing set of 2,032 pairs.

##Reference
> Han Yang, Ruimao Zhang, Xiaobao Guo, Wei Liu, and Wang
Zuo. Towards photo-realistic virtual try-on by adaptively generating-preserving image content. In IEEE Conf. on Computer Vision and Pattern Recognition (CVPR), pages 7850–7859. IEEE, 2020.
> Minqi Li, Richard Yida Xu, Jing Xin, Kaibing Zhang, and Junfeng Jing. Fast non-rigid points registration with cluster correspondences projection. Signal Processing, 170:107425, 2020