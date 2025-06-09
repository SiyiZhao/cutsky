/*******************************************************************************
* cutsky.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "convert_z.h"
#include "survey_geom.h"
#include "proc_cat.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
  /* Load configuration. */
  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return CUTSKY_ERR_CFG;
  }

  /* Prepare for redshift conversion. */
  ZCVT *zcvt;
  if (!(zcvt = zcvt_init(conf))) {
    printf(FMT_FAIL);
    P_EXT("failed to setup distance to redshift conversion\n");
    conf_destroy(conf);
    return CUTSKY_ERR_ZCVT;
  }

  /* Prepare for survey geometry application. */
  GEOM *geom;
  if (!(geom = geom_init(conf))) {
    printf(FMT_FAIL);
    P_EXT("failed to setup survey geometry\n");
    conf_destroy(conf);
    zcvt_destroy(zcvt);
    return CUTSKY_ERR_GEOM;
  }

  /* Run the cut-sky procedure. */
  if (process(conf, zcvt, geom)) {
    printf(FMT_FAIL);
    P_EXT("failed to generate the cut-sky catalog\n");
    conf_destroy(conf);
    zcvt_destroy(zcvt);
    geom_destroy(geom);
    return CUTSKY_ERR_CUTSKY;
  }

  conf_destroy(conf);
  zcvt_destroy(zcvt);
  geom_destroy(geom);
  return 0;
}
