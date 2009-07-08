
      PARAMETER       (MEGA=100 000 000) ! 100 mio words
*                           -----------
      COMMON/COMEGA/MQ(MEGA)
      REAL          QM(MEGA)
      PARAMETER          (MEGAH=MEGA/2)
      DOUBLE PRECISION DQ(MEGAH)
      EQUIVALENCE  (MQ(1),QM(1),DQ(1))

      COMMON/CMINDC/MQU,INDICS(0:39)
      EQUIVALENCE (INDICS( 0),IND0),(INDICS( 1),IND1),(INDICS( 2),IND2),
     +            (INDICS( 3),IND3),(INDICS( 4),IND4),(INDICS( 5),IND5),
     +            (INDICS( 6),IND6),(INDICS( 7),IND7),(INDICS( 8),IND8),
     +            (INDICS( 9),IND9)
      EQUIVALENCE (INDICS(10),INDA),(INDICS(11),INDB),(INDICS(12),INDC),
     +            (INDICS(13),INDD),(INDICS(14),INDE),(INDICS(15),INDF),
     +            (INDICS(16),INDG),(INDICS(17),INDH),(INDICS(18),INDI),
     +            (INDICS(19),INDJ),(INDICS(20),INDK),(INDICS(21),INDL),
     +            (INDICS(22),INDM),(INDICS(23),INDN),(INDICS(24),INDO),
     +            (INDICS(25),INDP),(INDICS(26),INDQ),(INDICS(27),INDR),
     +            (INDICS(28),INDS),(INDICS(29),INDT),(INDICS(30),INDU),
     +            (INDICS(31),INDV),(INDICS(32),INDW),(INDICS(33),INDX),
     +            (INDICS(34),INDY),(INDICS(35),INDZ)
      EQUIVALENCE (INDICS(36),IDOT1),(INDICS(37),IDOT2),
     +            (INDICS(38),IDOT3),(INDICS(39),IDOT4)

