# ゲノムのマッピング

- BWAを用いてリードを参照配列にマッピング
```
./bwa bwasw genomes/human shoga/FAH64410_caf194bd39f018b128560eba1ade5b8c7f6afc21_* > human+shoga.sam
```
- BAM形式に変換
```
samtools view -bS human+shoga.sam > human+shoga.bam
```
- 染色体位置でソートする
```
samtools sort human+shoga.bam human+shoga.sorted
```
- Indexを作成
```
samtools index human+shoga.sorted.bam
```
- マッピングの状況を見てみる
```
samtools flagstat human+shoga.sorted.bam
```
