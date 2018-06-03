__rethrow_casa_exceptions = True
h_init()
try:
    hif_restoredata (vis=['uid___A002_Xbbe66a_X1954'], session=['session_1'], ocorr_mode='ca')
finally:
    h_save()
