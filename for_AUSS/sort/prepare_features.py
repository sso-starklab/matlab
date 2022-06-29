import numpy as np
from scipy.spatial import distance
import time

# i,j - indexes
# next - reshape(make sure it works right), write the functions below

def prepare_features(i, j, clu, mean_spk, std_spk, cc, u_clu):
    s2 = time.time()
    idx1 = np.where(clu == u_clu[i])
    idx2 = np.where(clu == u_clu[j])

    s3 = time.time()
    mean_spk1 = mean_spk[:, :, i]
    [mean_spk1, ind] = trim_spk_8ch(mean_spk1)
    mean_spk1 = mean_spk1.flatten()
    mean_spk2 = (mean_spk[ind, :, j]).flatten()

    s4 = time.time()
    std_spk1 = (std_spk[ind, :, i]).flatten()
    std_spk2 = (std_spk[ind, :, j]).flatten()

    if np.sum(cc[:, i, i]) == 0:
        acc1 = cc[:, i, i]
    else:
        acc1 = cc[:, i, i] / np.max(cc[:, i, i])
    acc1 = acc1[20:]

    if np.sum(cc[:, j, j]) == 0:
        acc2 = cc[:, j, j]
    else:
        acc2 = cc[:, j, j] / np.max(cc[:, j, j])
    acc2 = acc1[20:]

    if np.sum(cc[:, i, j]) == 0:
        ccgtag = cc[:, i, j]
    else:
        ccgtag = cc[:, i, j] / np.max(cc[:, i, j])

    s8 = time.time()
    i1 = len(idx1[0])
    i2 = len(idx2[0])
    n = np.minimum(i1, i2) / np.maximum(i1, i2)

    s9 = time.time()

    x = x2feat1(mean_spk1, mean_spk2, std_spk1, std_spk2, acc1.T, acc2.T, ccgtag.T, n)
    x = featForAB2(x)

    # print( "Time 2: ", s3-s2, "Time 3: ", s4-s3, "Time 4: ", s5-s4, "Time 5: ", s6-s5, "Time 6: ", s7-s6, "Time 7: ",
    # s8-s7, "Time 8: ", s9-s8)

    return x


def trim_spk_8ch(mean_spk):
    n_channels = np.size(mean_spk, 0)
    max_idx = np.argmax(mean_spk[:, 15], 0)
    if max_idx <= 3:
        channels_idx = np.arange(0, 8)
    elif n_channels < max_idx + 5:
        channels_idx = np.arange(n_channels - 8, n_channels)
    else:
        channels_idx = np.arange(max_idx - 3, max_idx + 5)
    t8ch_mspk = mean_spk[channels_idx.T, :]
    return t8ch_mspk, channels_idx


def x2feat1(mean_spk1, mean_spk2, std_spk1, std_spk2, acc1, acc2, ccgtag, n):
    n_samp = 8 * 32

    # rearranging mean wave form and std
    all_spk = np.concatenate((mean_spk1, mean_spk2, std_spk1, std_spk2))
    mat_m_spk1 = np.reshape(mean_spk1, (8, 32)).T
    r = np.ptp(mat_m_spk1, axis=0)
    ch_order = np.argsort(r)
    ch_order = np.flip(ch_order)
    new_features = rearrange(all_spk, ch_order)

    # normalizing mean wave form and std
    m = np.max(abs(new_features[0:n_samp * 2]))
    new_features = new_features / m
    new_features = np.concatenate((new_features.T, acc1, acc2, ccgtag, np.array([n])))
    return new_features


def rearrange(all_spk, ch_order):
    mat = np.zeros((4, 1))
    for k in range(0, 8):
        s = 32 * (ch_order[k])
        e = s + 32
        c = 32 * 8
        wf1 = all_spk[s: e]
        wf2 = all_spk[s + c: e + c]
        sd1 = all_spk[s + (c * 2): e + (c * 2)]
        sd2 = all_spk[s + (c * 3): e + (c * 3)]

        mat2 = np.stack((wf1, wf2, sd1, sd2))
        mat = np.concatenate((mat, mat2), axis=1)
    mat = np.delete(mat, 0, axis=1)
    mat = mat.flatten()
    return mat


def featForAB2(x):
    dist = feat1(x)
    outSamp = feat2(x)
    Irclust = feat3(x)
    ccSamp = feat4(x)
    ccDist = feat5(x)
    peakDiff = feat6(x)
    nspkRation = feat7(x)
    spkSD = feat8(x)
    x_all_feat = np.concatenate((dist, outSamp, Irclust, ccSamp, ccDist, peakDiff, nspkRation, spkSD))
    return x_all_feat


# F1 euclidean dist between each channel
def feat1(x):
    dist1 = np.zeros(8)
    for j in range(0, 8):
        s = 32 * j
        e = s + 32
        c = 32 * 8

        wf1 = x[s: e]
        wf2 = x[s + c: e + c]
        dist1[j] = distance.euclidean(wf1, wf2)
    return dist1


# F2 number of samples outside of SD boundary per channel
def feat2(x):
    out_samp = np.zeros(8)
    for j in range(0, 8):
        s = 32 * j
        e = s + 32
        c = 32 * 8

        wf1 = x[s: e]
        wf2 = x[s + c: e + c]
        sd1 = x[s + (c * 2): e + (c * 2)]
        y = np.logical_or(wf2 > (wf1 + sd1), wf2 < (wf1 - sd1))
        y = np.count_nonzero(y) / 32
        out_samp[j] = y
    return out_samp


# F3 multiplication with max wf channel of the first cluster (inspired by ironclus)
def feat3(x):
    Irclust = np.zeros(8)
    mwf = x[0:32]
    for j in range(0, 8):
        s = 32 * j
        e = s + 32
        c = 32 * 8

        wf1 = x[s: e]
        wf2 = x[s + c: e + c]
        m1 = wf1 @ mwf.T
        m2 = wf2 @ mwf.T
        Irclust[j] = m1 - m2
    return Irclust


# F4 9 middle samples in the cch
def feat4(x):
    z = len(x) - 22
    cc_samp = x[z - 4:z + 5]
    return cc_samp


# F5 euclidean distance between acc1: acc2, acc1: ccg, acc2: ccg
def feat5(x):
    c = 32 * 8 * 4

    d1 = distance.euclidean(x[c:c + 21], x[c + 21:c + 42])
    d2 = distance.euclidean(x[c:c + 21], x[-22:-1])
    d3 = distance.euclidean(x[c + 21:c + 42], x[-22:-1])

    ccDist = np.array([d1, d2, d3])
    return ccDist


# F6 peak diff in each channel
def feat6(x):
    wf1 = x[0:255]
    wf2 = x[256:511]
    idxPeak = np.arange(15, 256, 32)
    peakDiff = abs(wf1[idxPeak] - wf2[idxPeak]) / abs(wf1[idxPeak] + wf2[idxPeak])
    peakDiff = np.array([sum(peakDiff)])
    return peakDiff


# F7 N spikes ratio
def feat7(x):
    return np.array([x[-1]])


# F8 4 first channels with their SD
def feat8(x):
    c = 32 * 4
    wf1 = x[:c]
    wf2 = x[c * 2: c * 3]
    sd1 = x[c * 4: c * 5]
    spkSD = np.concatenate((wf1, wf2, sd1))
    return spkSD
