// Make code cells collapsible and add global hide/show code button
document.addEventListener('DOMContentLoaded', function() {
    const inputWrappers = document.querySelectorAll('.jp-Cell-inputWrapper');
    console.log('Found ' + inputWrappers.length + ' input wrappers');
    
    // Create and insert the hide/show code button at the top
    const button = document.createElement('button');
    button.id = 'toggle-code-button';
    button.textContent = 'Hide Code';
    button.style.cssText = `
        position: fixed;
        top: 20px;
        left: 20px;
        z-index: 1000;
        padding: 10px 20px;
        background-color: #808080;
        color: white;
        border: none;
        border-radius: 5px;
        cursor: pointer;
        font-size: 14px;
        font-weight: 600;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        transition: background-color 0.2s ease;
    `;
    
    button.addEventListener('mouseover', function() {
        button.style.backgroundColor = '#606060';
    });
    
    button.addEventListener('mouseout', function() {
        button.style.backgroundColor = '#808080';
    });
    
    document.body.appendChild(button);
    
    let allCollapsed = false;
    
    button.addEventListener('click', function() {
        allCollapsed = !allCollapsed;
        button.textContent = allCollapsed ? 'Show Code' : 'Hide Code';
        
        inputWrappers.forEach((wrapper) => {
            if (!wrapper.parentElement.classList.contains('jp-CodeCell')) {
                return;
            }
            
            if (allCollapsed) {
                wrapper.classList.add('collapsed');
                wrapper.classList.remove('expanded');
            } else {
                wrapper.classList.remove('collapsed');
                wrapper.classList.add('expanded');
            }
        });
    });
    
    inputWrappers.forEach((wrapper, index) => {
        const codeCell = wrapper.parentElement;
        
        // Only handle code cells
        if (!codeCell.classList.contains('jp-CodeCell')) {
            return;
        }
        
        const editor = wrapper.querySelector('.jp-CodeMirrorEditor');
        if (!editor) return;
        
        // Start in expanded state
        wrapper.classList.remove('collapsed');
        wrapper.classList.add('expanded');
        
        // Add click handler
        wrapper.addEventListener('click', function(e) {
            // Ignore clicks on the code editor itself
            if (e.target.closest('.jp-CodeMirrorEditor')) {
                return;
            }
            
            e.preventDefault();
            e.stopPropagation();
            
            const isCollapsed = wrapper.classList.contains('collapsed');
            
            if (isCollapsed) {
                wrapper.classList.remove('collapsed');
                wrapper.classList.add('expanded');
            } else {
                wrapper.classList.add('collapsed');
                wrapper.classList.remove('expanded');
            }
            
            return false;
        }, true);
    });
});

