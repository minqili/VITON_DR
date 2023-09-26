import cv2
import numpy as np
import scipy.io as sio
from scipy.spatial.distance import pdist
def curve(x2, y2):
    curveture = []
    for i in range(len(x2)-2):
        x=(x2[i]-x2[i-1],y2[i]-y2[i-1])
        y=(x2[i+1]-x2[i],y2[i+1]-y2[i])
        d = 1 - pdist([x, y], 'cosine')
        sin = np.sqrt(1-d**2)
        dis = np.sqrt((x2[i-1]-x2[i+1])**2 + (y2[i-1]-y2[i+1])**2)
        k = 2 * sin / dis
        curveture.append(k[0])
    return curveture

def get_cloth_data(Image):
    blurred = cv2.GaussianBlur(Image, (7, 7), 0)
    gaussImg = cv2.Canny(blurred, 20, 20)
    W, H = gaussImg.shape
    list_i = []
    list_j = []
    p = 0
    p_1 = 0
    img_zeros = np.zeros((256, 192), np.uint8)
    for i in range(W):
        for j in range(H):
            if gaussImg[i][j] == 255 and p % 1 == 0:
                img_zeros[i][j] = gaussImg[i][j]
                list_i.append((-i + 255) / 255)
                list_j.append((j) / 191)
                p_1 += 1  
            p += 1  
    Keypoint = np.array(list(zip(list_j, list_i)))
    Curve = curve(Keypoint[:, 0], Keypoint[:, 1])
    KP=[]
    for s in range(0, int((Keypoint.shape[0])/15)*15, 15):
        tem_curve = Curve[s:s+3]
        max_v = max(tem_curve)
        id_max = tem_curve.index(max_v)
        k_t = Keypoint[s:s+3]
        KP.append(k_t[id_max])
    KP =np.array(KP)
    file_name = ".\\models\\datamat\\data1.mat"
    data = sio.loadmat(file_name)
    data["x1"] = KP
    sio.savemat(".\\models\\datamat\\data1.mat", data)

def get_mask_data(Image):
    blurred = cv2.GaussianBlur(Image, (7, 7), 0)
    gaussImg = cv2.Canny(blurred, 20, 20)
    W, H = gaussImg.shape
    list_i = []
    list_j = []
    p = 0
    p_1 = 0
    img_zeros = np.zeros((256, 192), np.uint8)
    for i in range(W):
        for j in range(H):
            if gaussImg[i][j] == 255 and p % 1 == 0:
                img_zeros[i][j] = gaussImg[i][j]
                list_i.append((-i + 255) / 255)
                list_j.append((j) / 191)
                p_1 += 1
            p += 1
    Keypoint = np.array(list(zip(list_j, list_i)))
    Curve = curve(Keypoint[:, 0], Keypoint[:, 1])
    KP=[]
    for s in range(0, int((Keypoint.shape[0])/15)*15, 15):
        tem_curve = Curve[s:s+3]
        max_v = max(tem_curve)
        id_max = tem_curve.index(max_v)
        k_t = Keypoint[s:s+3]
        KP.append(k_t[id_max])
    KP =np.array(KP)
    file_name = ".\\models\\datamat\\data1.mat"
    data = sio.loadmat(file_name)
    data["y2a"] = KP
    sio.savemat(".\\models\\datamat\\data1.mat", data)