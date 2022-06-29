% cell array of sessions from CA1
regCA1      = { 'mS234_19', 'mA234_15',  'mC41_31', 'mC41_32', 'mC41_33', 'mC41_34', ...
    'mC41_38', 'mC41_39', 'mC41_40', 'mC41_41', 'mC41_43', 'mC41_44', 'mC41_45', 'mC41_53', 'mC41_55',...
   'mP101_06', 'mP101_07', 'mP101_10', 'mP101_12', 'mP101_13', 'mP101_14', 'mP23_12', 'mP23_13', 'mP23_14', 'mP23_15',...
  'mP23_16', 'mP23_18', 'mP23_20', 'mP23_22', 'mP23_26', 'mP23_27', 'mP23_29', 'mP23_32', 'mP23_33', 'mP23_34', 'mP23_35' ...
'mP23_36', 'mP23_37', 'mP23_38', 'mP23_39', 'mP23_41', 'mP23_42' };
% cell array of sessions from nCX
nCX         = { 'mA234_01', 'mA234_14', 'mA234_15', 'mA234_16', 'mB142_06', 'mC41_10', 'mC41_11', 'mC41_12', 'mC41_13',...
    'mDL5_05', 'mDL5_06', 'mF105_10', 'mF108_01', 'mF108_02', 'mF79_26', 'mF84_01', 'mF84_02', 'mF84_03', 'mF84_04', ...
    'mF84_05', 'mF84_06', 'mF84_08', 'mF84_09', 'mF93_02', 'mF93_04', 'mP23_03', 'mP23_04', 'mS234_01', 'mS234_02', 'mV99_03'};
% cell array of sessions from HIP
HIP         = {'mB142_09', 'mC41_19', 'mC41_20'};
% cell array of sessions from WM
WM          = { 'mP23_05', 'mP23_06'};
% cell array of sessions from nCX layers 4-5 + WM + CA1
mixed       = { 'mV99_07', 'mV99_09', 'mV99_10', 'mK01_10', 'mK01_13', 'mK01_14'};
% cell array of sessions from CA1 + lower areas
CA1_DG      = { 'mV99_11', 'mV99_12', 'mV99_13', 'mV99_14', 'mV99_15', 'mV99_16', 'mV99_17', 'mV99_18', 'mV99_19', 'mK01_16' };


% cell array of mouse lines
VIP         = {'mV99'};
CamK        = {'mF84', 'mF79', 'mF93', 'mF108'};
PV          = {'mP23', 'mP101', 'mDL5'};
FVB_CamK    = {'mC41', 'mA234', 'mS234'};
CCK         = {'mK01'};
FVB_c57B    = {'mF105'};
FVB_PV      = {'mB142'};
