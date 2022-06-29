from prepare_features import *

clu = np.load('/home/tali/matlab/AUSS_python/mP31_04.clu.1.npy')
mspk = np.load('/home/tali/matlab/AUSS_python/mP31_04.mspk.1.npy')
sspk = np.load('/home/tali/matlab/AUSS_python/mP31_04.sspk.1.npy')
nspk_vec = np.load('/home/tali/matlab/AUSS_python/mP31_04.nspk_vec.1.npy')
cc = np.load('/home/tali/matlab/AUSS_python/mP31_04.cc.1.npy')
nspk_vec = np.squeeze(nspk_vec)
#f1 = np.load('/home/tali/matlab/AUSS_python/mP31_04/mP31_04.X2Feat1.1-1-1.npy')
f2 =np.load('/home/tali/matlab/AUSS_python/mP31_04/mP31_04.featForAB2.1-1-2.npy').squeeze()
i = 0
j = 2
x = prepare_features(i, j, clu, mspk, sspk, cc)
bool = (x==f2)
print(bool)
