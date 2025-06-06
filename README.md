# LAMMPS Bulk Job Sample

## 概要

物性研究所システムBにおいてLAMMPSをバルクジョブで実行するためのサンプルファイル。

テストのためにi8cpuキューを使っているが、実際のプロダクトランはFキューなどを使うこと。

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
    cd data
    sbatch job1.sh
    sbatch job8.sh
    ```

## ジョブの内容

50万原子を1ノード128並列(flat-MPI)で10000ステップ実行するLAMMPSのジョブである。T=0.90からT=1.04まで、異なる温度で温度制御をしている。

シリアルなら

```sh
srun lammps < T090/T090.input
```

などとして実行するジョブである。それぞれ概ね140秒程度のジョブであり、ファイルの読み込みなども含めても200秒程度で終わるはずである。

実行終了後、ログ(`slurm-???????.out`)を確認すると、それぞれの実行時間がわかる。


シリアルジョブなら

```txt
Elapsed time: 1374 seconds
```

と、全体で1374秒かかっているのに対して、

例えば、バルクジョブなら

```txt
Elapsed time: 213 seconds
```

と、大幅に実行時間が短縮していることがわかる。

それぞれのディレクトリに`*.log`ファイルが生成される。例えば`T090/T090.log`なら

```txt
   Step          Temp          E_pair         E_mol          TotEng         Press     
         0   0.33335784    -8.1018861      0             -7.6018503     -0.038231056  
      1000   0.35872634    -7.5786182      0             -7.0405298      1.5181176    
      2000   0.81253608    -6.9566444      0             -5.7378427      2.4099187    
      3000   0.90345209    -6.641065       0             -5.2858895      0.76827878   
      4000   0.90169206    -6.499719       0             -5.1471836      0.39242165   
      5000   0.90010336    -6.4356147      0             -5.0854624      0.46632326   
      6000   0.89993486    -6.3759544      0             -5.0260548      0.75303208   
      7000   0.90065031    -6.2347114      0             -4.8837386      1.7027286    
      8000   0.90110407    -6.1486899      0             -4.7970365      2.4674989    
      9000   0.8975787     -6.2253512      0             -4.8789858      2.1938972    
     10000   0.8986377     -6.298002       0             -4.9500482      1.6930766   
```

と、温度がおよそ0.9に制御されており、`T092/T092.log`なら

```txt
   Step          Temp          E_pair         E_mol          TotEng         Press     
         0   0.33335784    -8.1018861      0             -7.6018503     -0.038231056  
      1000   0.36006527    -7.576788       0             -7.0366912      1.5255005    
      2000   0.82481222    -6.9431712      0             -5.7059553      2.4547777    
      3000   0.92389541    -6.615707       0             -5.2298667      0.83335939   
      4000   0.92068515    -6.4710826      0             -5.0900576      0.52411015   
      5000   0.91890035    -6.4081691      0             -5.0298213      0.61943105   
      6000   0.91906177    -6.3480456      0             -4.9694557      0.92308146   
      7000   0.92076604    -6.1914581      0             -4.8103118      1.9642132    
      8000   0.92101113    -6.11327        0             -4.731756       2.6787946    
      9000   0.92255472    -6.1957347      0             -4.8119054      2.3633792    
     10000   0.92107821    -6.2608367      0             -4.8792222      1.8851376    
```

と、温度が0.92に制御されていることがわかる。

## ディレクトリ構造

`prepare.py`を実行すると、以下のようなディレクトリ構造となる。

```txt
├── data
│   ├── job1.sh
│   ├── job8.sh
│   ├── T090
│   │   ├── input.atoms
│   │   └── T090.input
│   ├── T092
│   │   ├── input.atoms
│   │   └── T092.input
│   ├── T094
│   │   ├── input.atoms
│   │   └── T094.input
│   ├── T096
│   │   ├── input.atoms
│   │   └── T096.input
│   ├── T098
│   │   ├── input.atoms
│   │   └── T098.input
│   ├── T100
│   │   ├── input.atoms
│   │   └── T100.input
│   ├── T102
│   │   ├── input.atoms
│   │   └── T102.input
│   └── T104
│       ├── input.atoms
│       └── T104.input
```

すなわち、`data`ディレクトリ直下にジョブスクリプトが、その下にパラメータごとにディレクトリが作成され、それぞれのディレクトリに初期位置ファイル(`*.atoms`)、インプットファイル(`*.input`)が配置される。

例えば`T090`には温度T=0.90のジョブが配置されている。しかし、ジョブの実行をその一つ上のディレクトリ(`data`ディレクトリ)で実行する関係から、原子ファイルとログファイルの位置がずれることに注意。具体的には`T090.input`は以下のような内容となる。

```txt
units lj
atom_style atomic
boundary p p p
timestep 0.001

read_data T090/input.atoms

mass 1 1.0

pair_style lj/cut 3.0
pair_coeff 1 1 1.0 1.0 3.0

fix 1 all nvt temp 0.9 0.9 0.1
log T090/T090.log
thermo 1000
run 10000
```

`read_data`や`log`コマンドが、`T090/`というディレクトリを指示している。これはジョブが一つ上のディレクトリから実行されるから(カレントディレクトリが一つ上にずれるから)である。

`prepare.py`は、それぞれの温度で`generate_config.py`を実行するラッパースクリプトであり、実際のファイルは`generate_config.py`が生成している。`generate_config.py`は、密度と温度を指定すると、温度に対応するディレクトリを作成して原子ファイルとインプットファイルを作成するスクリプトであり、例えば、

```sh
python3 generate_config.py -d 0.9 -t 0.9
```

を実行すると、`data`ディレクトリ以下に(もし存在しなければ)`T090`ディレクトリを作成し、その下に

```txt
data/T090/input.atoms
data/T090/T090.input
```

を生成する。今回は簡単のために密度は全て0.9としたが、密度も変更するならディレクトリ名に含めると良い。

`prepare.py`はジョブスクリプトも生成する。`job1.sh`がシリアルジョブであり、`job8.sh`がバルクジョブである。

シリアルジョブは以下のようなファイルとなる。

```sh:job1.sh
#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 128

source /home/issp/materiapps/intel/lammps/lammpsvars.sh
SECONDS=0

srun lammps < T090/T090.input
srun lammps < T092/T092.input
srun lammps < T094/T094.input
srun lammps < T096/T096.input
srun lammps < T098/T098.input
srun lammps < T100/T100.input
srun lammps < T102/T102.input
srun lammps < T104/T104.input

echo "Elapsed time: ${SECONDS} seconds"
```

`SECONDS`を用いて実行時間を測定している。

バルクジョブ用のジョブスクリプト`job8.sh`は以下の通り。

```sh:job8.sh
#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 8


source /home/issp/materiapps/intel/lammps/lammpsvars.sh
SECONDS=0

srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T090/T090.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T092/T092.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T094/T094.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T096/T096.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T098/T098.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T100/T100.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T102/T102.input &
sleep 5
srun --exclusive --mem-per-cpu=1840 -N 1 -n 128 -c 1 lammps < T104/T104.input &
sleep 5

wait
echo "Elapsed time: ${SECONDS} seconds"
```

`srun`に、`--exclusive`を指定して計算資源の排他的占有を指示する。また、`-N 1 -n 128 -c 1`により、1ノード、128プロセス、1スレッド/1プロセス、すなわち128プロセスのFlat-MPIジョブとして実行している。`sleep 5`コマンドを挟んでいるのは、経験上その方が安定するからである(たまに計算資源の排他占有に失敗する？)。やはり最後に`SECONDS`を利用して実行時間を測定している。

## ライセンス

MIT
