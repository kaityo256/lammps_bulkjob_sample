# LAMMPS Bulk Job Sample

## 概要

物性研究所システムBにおいてLAMMPSをバルクジョブで実行するためのサンプルファイルです。

テストのためにi8cpuキューを使っていますが、実際のプロダクトランはFキューなどを使ってください。

## 使い方

1. リポジトリをクローンする。
    ```sh
    git clone git@github.com:kaityo256/lammps_bulkjob_sample.git
    cd lammps_bulkjob_sample
    ```
1. ジョブスクリプトやインプットファイル、原子ファイルを生成する。
    ```sh
    python3 prepare.py
    ```
1. `data`ディレクトリに、8本のジョブをシリアルに実行する`job1.sh`と、バルクジョブとしてまとめて実行する`job8.sh`ができているので、ジョブを投入する。
    ```sh
    sbatch job1.sh
    sbatch job8.sh
    ```

実行終了後、ログを確認すると、それぞれの実行時間がわかる。

## ライセンス

MIT
