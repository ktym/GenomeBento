# ゲノムのマッピング

## BWAを用いてリードを参照配列にマッピング

```sh
./bwa bwasw databases/genomes/human samples/Shoga/gDNA/FAH64410_caf194bd39f018b128560eba1ade5b8c7f6afc21_* > mappings/human+shoga/human+shoga_gDNA.sam
```

### ここでディレクトリを移動

```sh
cd appings/human+shoga/
```

## BAM形式に変換

```sh
samtools view -bS human+shoga_gDNA.sam > human+shoga_gDNA.bam
```

## 染色体位置でソートする

```sh
samtools sort human+shoga_gDNA.bam human+shoga_gDNA.sorted
```
## Indexを作成

```sh
samtools index human+shoga_gDNA.sorted.bam
```

## マッピングの状況を見てみる

```sh
samtools flagstat human+shoga_gDNA.sorted.bam
```
