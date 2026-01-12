#!/usr/bin/env python3
"""
Gemini API 测试脚本
用于验证 Gemini API 配置是否正确
"""

import os
import google.generativeai as genai
from pathlib import Path

def load_api_key():
    """从环境变量或 .env 文件加载 API Key"""
    # 首先尝试从环境变量获取
    api_key = os.getenv('GEMINI_API_KEY')

    if not api_key:
        # 尝试从 .env 文件读取
        env_file = Path(__file__).parent / '.env'
        if env_file.exists():
            with open(env_file, 'r') as f:
                for line in f:
                    if line.startswith('GEMINI_API_KEY='):
                        api_key = line.strip().split('=', 1)[1]
                        break

    return api_key

def test_gemini_api():
    """测试 Gemini API 连接"""
    print("=" * 50)
    print("Gemini API 配置测试")
    print("=" * 50)

    # 加载 API Key
    api_key = load_api_key()

    if not api_key or api_key == 'your_api_key_here':
        print("\n❌ 错误: 未找到有效的 API Key")
        print("\n请按以下步骤配置:")
        print("1. 访问 https://makersuite.google.com/app/apikey 获取 API Key")
        print("2. 复制 .env.example 为 .env: cp .env.example .env")
        print("3. 编辑 .env 文件，将 your_api_key_here 替换为你的实际 API Key")
        print("4. 或者设置环境变量: export GEMINI_API_KEY='your_api_key'")
        return False

    try:
        # 配置 API
        genai.configure(api_key=api_key)

        # 列出可用模型
        print("\n✓ API Key 配置成功")
        print("\n可用的模型:")
        for model in genai.list_models():
            if 'generateContent' in model.supported_generation_methods:
                print(f"  - {model.name}")

        # 测试简单的文本生成
        print("\n正在测试文本生成...")
        model = genai.GenerativeModel('gemini-pro')
        response = model.generate_content("Say 'Hello, Gemini API is working!' in Chinese")

        print("\n✓ 测试成功!")
        print(f"\n模型响应: {response.text}")

        return True

    except Exception as e:
        print(f"\n❌ 错误: {str(e)}")
        print("\n请检查:")
        print("1. API Key 是否正确")
        print("2. 网络连接是否正常")
        print("3. API Key 是否有足够的配额")
        return False

if __name__ == "__main__":
    test_gemini_api()
