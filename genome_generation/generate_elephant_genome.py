import random
def generate_artificial_elephant_genome(length=3000000000):
    """
    코끼리 유전체처럼 보이는 염기서열을 인위적으로 생성한다.
    Parameters:
        length (int): 생성할 유전체 길이
    Returns:
        str: 랜덤으로 생성된 ACGT 시퀀스 문자열
    """
    bases = ['A', 'C', 'G', 'T']
    genome = ''.join(random.choices(bases, k=length))
    return genome

def save_reference_as_txt(sequence, filename="reference_100K.txt"):
    with open(filename, "w") as f:
        f.write(sequence)


# 100만 bp짜리 유전체 생성
genome_seq = generate_artificial_elephant_genome(length=100000)

save_reference_as_txt(genome_seq)