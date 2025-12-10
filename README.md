# exCGen_MCP

An MCP server implementation of exCGen using FastMCP. Ported from the original web version: https://github.com/kkomdori/exCGen_v1
  
### Features of exCGen_MCP:
**Design mutually exclusive DNA Sequences:**
- Generate mutually exclusive DNA sequences that won't bind to each other
- Specify custom length (in base pairs)
- Control GC content (percentage of guanine and cytosine)
- Create multiple sequence pairs at once
- Use preset patterns with specific sequences embedded
- Avoid specific sequences (like restriction sites)
- Ensure new sequences don't bind to existing ones

**Use Cases:**
- Designing primers for PCR
- Creating molecular barcodes
- Developing DNA probes
- Synthetic biology applications

# Connecting this MCP to Claude Desktop in local condition

This guide explains how to integrate the **exCGen MCP** server with the Claude Desktop App.   
This project uses **`uv`** (an extremely fast Python package manager) to handle dependencies and execution automatically.

### 1. Prerequisites

Before you begin, ensure you have the following installed:

1. **Claude Desktop App**: [Download here](https://claude.ai/download). (Tested Version: 1.0.1768)
2. **Git**: [Download here](https://git-scm.com/) (Required to download the code).
3. **uv (Required)**: This tool manages the Python environment for this project.
    
    - **Windows (PowerShell)**:
        
        PowerShell
        ```
        powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
        ```
        
    - **macOS / Linux**:
        
        Bash
        ```
        curl -LsSf https://astral.sh/uv/install.sh | sh
        ```
      
    > **Note:** After installation, restart your terminal and type `uv --version` to verify it is installed correctly.  

### 2. Installation

Open your terminal (Command Prompt or PowerShell), navigate to the folder where you want to save the project, and run the following commands:

Bash
```
git clone https://github.com/kkomdori/exCGen_MCP.git
cd exCGen_MCP
```

_(You do not need to manually install dependencies. `uv` will handle them automatically when connected.)_

### 3. Locate Configuration File

You need to edit the Claude Desktop configuration file to tell it where to find this project.

1. **Open the configuration file location**:
    - **Windows**: Press `Win + R`, paste the following, and press Enter: `%APPDATA%\Claude\`
    - **macOS**: Open Finder, press `Cmd + Shift + G`, and enter: `~/Library/Application Support/Claude/`
        
2. Look for a file named `claude_desktop_config.json`.
    - If it doesn't exist, create a new file with that exact name.
    - Open it with a text editor (Notepad, VS Code, etc.).
        
### 4. Edit Configuration (Crucial Step)

Copy and paste the code below into your `claude_desktop_config.json` file. **⚠️ IMPORTANT:** You must replace the placeholders (e.g., `[YOUR_USERNAME]`) with your actual Windows username.
This configuration uses `uv` to run the project and includes specific environment variables to prevent encoding issues (text display errors).

JSON
```
{
  "mcpServers": {
    "exCGen-MCP": {
      "command": "C:\\Users\\[YOUR_USERNAME]\\.local\\bin\\uv.exe",
      "args": [
        "--directory",
        "C:\\Users\\[YOUR_USERNAME]\\Github\\exCGen_MCP",
        "run",
        "main.py"
      ],
      "env": {
        "PYTHONUTF8": "1",
        "PYTHONIOENCODING": "utf-8"
      }
    }
  }
}
```

### What you need to change:

1. **`command`**: The path to your `uv.exe`.
    - Typically: `C:\\Users\\YOUR_NAME\\.local\\bin\\uv.exe`
        
2. **`args` (Project Path)**: The path to where you downloaded this project.
    - Verify that `main.py` exists inside the folder you specify.
        
3. **JSON Formatting**:
    - On Windows, you must use **double backslashes** (`\\`) for paths.
    - Example: `C:\\Users\\JohnDoe\\...`

### 5. Verify Connection

1. **Restart** the Claude Desktop App completely.
2. Look for the **Plug icon (🔌)** on the right side of the input bar.
3. Click it and check if `exCGen-MCP` is listed with a green **Connected** status.
