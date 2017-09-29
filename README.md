cornet
================

インストール
------------

``` r
#install.packages("devtools") # もしなければ
devtools::install_github("shkonishi/cornet")
```

    ## Skipping install of 'cornet' from a github remote, the SHA1 (16ffbf21) has not changed since last install.
    ##   Use `force = TRUE` to force installation

関数及び外部データ一覧
----------------------

``` r
library(cornet)
ls("package:cornet")
```

    ## [1] "cluster_dat"  "cluster_mat"  "cluster_mine" "corheat"     
    ## [5] "matoedge"

### data

``` r
# data: normalized fpkm
fp <- system.file("extdata/nfpkm_200.txt", package = "cornet")
dat <- read.table(fp, header=TRUE, stringsAsFactors = FALSE)

# 200 genes
dat[1:6,1:6]; dim(dat)
```

    ##   id runs days reps  gene266  gene372
    ## 1  1    1    3    1 17.18686 4.130407
    ## 2  2    1    3    1 19.13915 4.413801
    ## 3  3    1    3    1 19.14471 3.961325
    ## 4  4    1    4    2 15.22853 4.668531
    ## 5  5    1    4    2 18.55824 6.054861
    ## 6  6    1    4    2 18.79719 7.802054

    ## [1] 108 204

### cluster\_mat

-   `amap::Dist`のメソッドから距離定義を選択
-   別手法で作成した距離行列を`as.dist`で変換したdistオブジェクトでも良い

``` r
res <- cornet::cluster_mat(dat = dat, distm = "spearman", clm = "average",
                           column = 5:ncol(dat), method_dycut = "hybrid",
                           x_fctr = dat$days, y_fctr = dat$runs, rep_fctr = dat$reps)
```

    ##  ..cutHeight not given, setting it to 268000  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-1.png)

``` r
# cutreeDynamicの結果
head(res$dynamic_cut)
```

    ## gene266 gene372 gene572 gene906 gene201 gene894 
    ##       1       2       2       2       1       1

``` r
# クラスタ別のデータフレーム
sapply(res$cluster_dat, dim)
```

    ##        1   2   3
    ## [1,] 108 108 108
    ## [2,]  96  82  22

``` r
# クラスタ別の
res$gg_mat_all
```

![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-2.png)

``` r
res$gg_mat_med
```

![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-3.png)

### cluster\_mine

-   `minerva::mine`の出力を整形, pearsonとspearmanも加える

``` r
# mineを連続実行、結果を整形出力
res.clm <- cluster_mine(cl_dat = res$cluster_dat)

# 
lapply(res.clm, function(x)x[1:6,])
```

    ## $`1`
    ##         x_id    y_id mic        mas mev      mcn      micr2      gmic
    ## 4425 gene204 gene655   1 0.03285807   1 2.000000 0.05284526 1.0000000
    ## 4485  gene67 gene519   1 0.01908699   1 3.000000 0.05036932 0.9935214
    ## 4460 gene179 gene519   1 0.04402342   1 3.000000 0.06798504 0.9942813
    ## 4470 gene789  gene67   1 0.04047976   1 3.584963 0.05791981 0.9806749
    ## 4457 gene179  gene67   1 0.01189086   1 3.000000 0.09046485 0.9930227
    ## 4484  gene67 gene371   1 0.09622638   1 3.000000 0.02788684 0.9865954
    ##           tic   pearson  spearman
    ## 4425 18.12556 0.9732188 0.9794221
    ## 4485 17.94552 0.9744900 0.9627502
    ## 4460 17.61939 0.9654092 0.9593205
    ## 4470 17.51823 0.9706082 0.9639506
    ## 4457 17.49004 0.9536955 0.9626930
    ## 4484 17.36597 0.9859580 0.9648080
    ## 
    ## $`2`
    ##         x_id    y_id       mic        mas       mev      mcn       micr2
    ## 3170 gene850 gene558 0.9744918 0.02346228 0.9744918 3.807355 0.008202131
    ## 3186 gene491 gene558 1.0000000 0.01791730 1.0000000 4.000000 0.059638626
    ## 3169 gene850 gene491 0.9744918 0.06201446 0.9744918 3.321928 0.052296549
    ## 2871 gene507 gene558 1.0000000 0.09008929 1.0000000 3.000000 0.250265011
    ## 2999 gene572 gene900 1.0000000 0.03026810 1.0000000 3.584963 0.141137013
    ## 3135 gene651 gene558 1.0000000 0.03673405 1.0000000 4.000000 0.123150137
    ##           gmic      tic   pearson  spearman
    ## 3170 0.9343362 17.28197 0.9830003 0.9613498
    ## 3186 0.9452102 16.93374 0.9697223 0.9770309
    ## 3169 0.9407437 16.76802 0.9603100 0.9584441
    ## 2871 0.9935214 16.48217 0.8658724 0.8620424
    ## 2999 0.9611421 16.45482 0.9267486 0.8917850
    ## 3135 0.9404114 16.26148 0.9364026 0.9300066
    ## 
    ## $`3`
    ##        x_id    y_id       mic        mas       mev      mcn      micr2
    ## 204 gene380 gene132 0.7558046 0.13125066 0.7558046 2.584963 0.06480325
    ## 112 gene530 gene498 0.6459352 0.11888400 0.6459352 2.000000 0.10039683
    ## 197 gene398 gene132 0.7145259 0.22761319 0.7145259 2.584963 0.19113752
    ## 189 gene917 gene132 0.5175286 0.07059585 0.5175286 2.000000 0.02622865
    ## 188 gene917 gene380 0.4971026 0.09883919 0.4971026 2.000000 0.04236234
    ## 180 gene805 gene132 0.5859004 0.15385737 0.5859004 2.584963 0.15300142
    ##          gmic       tic   pearson  spearman
    ## 204 0.6789852 11.021268 0.8312649 0.8333095
    ## 112 0.5718684  9.450314 0.7386057 0.7332495
    ## 197 0.6046873  9.192770 0.7234558 0.7442625
    ## 189 0.4661657  7.869230 0.7009279 0.6823287
    ## 188 0.4671553  7.828284 0.6743443 0.6715158
    ## 180 0.5013435  7.721901 0.6579506 0.6860061

### matoedge

``` r
# 相関行列のような対称行列から重み付きエッジリストを作成
edge.list <- matoedge(cor(res$cluster_dat[[1]]))
head(edge.list)
```

    ##      x_id    y_id       value
    ## 1 gene478 gene315 0.371799358
    ## 2 gene478 gene438 0.216129367
    ## 3 gene478 gene324 0.419083487
    ## 4 gene478 gene707 0.077285184
    ## 5 gene478  gene23 0.009250348
    ## 6 gene478 gene200 0.200931931
