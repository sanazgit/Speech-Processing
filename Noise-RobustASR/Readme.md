# Noise-Robust ASR with Whisper and Wav2Vec2

This project demonstrates a practical and modular pipeline for **automatic speech recognition (ASR)** in real-world, noisy conditions.  
It compares the performance of two leading ASR models — **OpenAI Whisper** and **Facebook Wav2Vec2** — under both raw and denoised audio inputs.

---

##  Features

-  Denoising raw audio using `noisereduce`
-  Visualization of audio features: Spectrogram and MFCC
-  Speech-to-text transcription using:
  - OpenAI Whisper
  - Facebook Wav2Vec2 (via HuggingFace)
-  Performance evaluation using **Word Error Rate (WER)**

---

## How to Run

1. Open the notebook in [Google Colab](https://colab.research.google.com/)
2. Install the required dependencies
3. Run each cell step-by-step
4. Observe outputs and comparisons

##  Result Summary

| Model     | Audio Type | Transcription                          |
|-----------|------------|----------------------------------------|
| Whisper   | Raw        | Why should one halt on the way?       |
| Whisper   | Denoised   | Wash it one halt on the way. ❌       |
| Wav2Vec2  | Raw        | WHY SHOULD ONE HALT ON THE WAY ✅     |
| Wav2Vec2  | Denoised   | WHY SHOULD ONE HALT ON THE WAY ✅     |

>  **Wav2Vec2 shows more robustness to noise compared to Whisper in this case.**


## Tech Stack

- Python 3.x
- PyTorch
- HuggingFace Transformers
- `openai-whisper`
- `torchaudio`, `librosa`, `noisereduce`


